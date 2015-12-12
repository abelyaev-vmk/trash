#include "ClothSim.h"
#include <cstdint>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


const float ClothMeshData::InitialHardness = 100;
const int ClothMeshData::elements = 34;
const float ClothMeshData::g = 0.6;
const float ClothMeshData::connectivity = 2.;
const float ClothMeshData::time_scale = 0.1;
const float ClothMeshData::wind = 2;
const float ClothMeshData::edge_length = 0.02;
const float ClothMeshData::loss = 0.80;

float dist(float4 a, float4 b) {
	return sqrt(dot((a - b),(a - b)));
}
ClothMeshData CreateTest2Vertices()
{
  ClothMeshData mdata;
  int el = ClothMeshData::elements;
  float edgeLength = ClothMeshData::edge_length;
  float hardnes = ClothMeshData::InitialHardness;
  mdata.vertPos0.resize(0);
  mdata.vertVel0.resize(0);
  mdata.vertPos1.resize(0);
  mdata.vertVel1.resize(0);
  mdata.texCoord.resize(0);
  mdata.vertForces.resize(0);

  mdata.edgeIndices.resize(0);
  mdata.edgeHardness.resize(0);
  mdata.edgeInitialLen.resize(0);
  mdata.vertNormals.resize(0);
  mdata.isFixed.resize(0);
  mdata.triangles.resize(0);

  for (int i = 0; i < el; ++i)
	  for (int j = 0; j < el; ++j) { 
		  mdata.vertPos0.push_back(float4(-0.5*el*edgeLength + i * edgeLength, -0.5*el*edgeLength + j * edgeLength, 0, 1));
		  mdata.vertVel0.push_back(float4(0, 0, 0, 0));
		  mdata.isFixed.push_back(j == (el - 1));
		  mdata.vertForces.push_back(float4(0, 0, 0, 0));
		  mdata.vertNormals.push_back(float3(0, 1, 0));
		  mdata.texCoord.push_back(float2((float)i/el, (float)j/el));
	  }

  for (int i = 0; i < el - 1; ++i)
	  for (int j = 0; j < el - 1; ++j) {
		  // left
		  mdata.triangles.push_back(i*el + j);
		  mdata.triangles.push_back(i*el + j + 1);
		  mdata.triangles.push_back((i + 1)*el + j + 1);

		  // rigth
		  mdata.triangles.push_back(i*el + j);
		  mdata.triangles.push_back((i + 1)*el + j + 1);
		  mdata.triangles.push_back((i + 1)*el + j);
	  }

  mdata.vertPos1 = mdata.vertPos0;
  mdata.vertVel1 = mdata.vertVel0;

  for (int i = 0; i < mdata.vertPos0.size(); ++i)
	  for (int j = i+1; j < mdata.vertPos0.size(); ++j) {
		  float d = dist(mdata.vertPos0[i], mdata.vertPos0[j]);
		  if (d <= ClothMeshData::connectivity*sqrt(2.0)*edgeLength) {
			  mdata.edgeIndices.push_back(i);
			  mdata.edgeIndices.push_back(j);
			  mdata.edgeHardness.push_back(hardnes);
			  mdata.edgeInitialLen.push_back(d);
		  }
	  }
  std::cout << "Total nodes:" << mdata.vertPos0.size() << std::endl;
  std::cout << "Total edges:" << mdata.edgeIndices.size() / 2 << std::endl;
  mdata.global_time = 0.0f;
  mdata.g_wind = float4(0, 0, 0, 0);

  mdata.texture = CreateGLTextureFromFile("../data/sponza_fabric_green_diff.tga");

  // you can use any intermediate mesh representation or load data to GPU (in VBOs) here immediately.                              <<===== !!!!!!!!!!!!!!!!!!

  // create graphics mesh; SimpleMesh uses GLUS Shape to store geometry; 
  // we copy data to GLUS Shape, and then these data will be copyed later from GLUS shape to GPU 
  //
  mdata.pMesh = std::make_shared<SimpleMesh>();


  // for tri mesh you will need normals, texCoords and different indices
  // 

  

  return mdata;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ClothMeshData::updatePositionsGPU(int type_)
{
  if (pMesh == nullptr)
    return;
 // copy current vertex positions to positions VBO

  GLUSshape& shape = pMesh->m_glusShape;
  shape.numberVertices = vertPos0.size();
  shape.vertices = (GLUSfloat*)malloc(4 * shape.numberVertices * sizeof(GLUSfloat));
  float4 *data = (pinPong ? &vertPos1[0] : &vertPos0[0]);
  memcpy(pMesh->m_glusShape.vertices, data, sizeof(float) * 4 * shape.numberVertices);
  shape.texCoords = (GLUSfloat*)malloc(2 * shape.numberVertices * sizeof(GLUSfloat));
  memcpy(pMesh->m_glusShape.texCoords, &texCoord[0], sizeof(float) * 2 * shape.numberVertices);

  if (type_ == 1) {
	  shape.numberIndices = edgeIndices.size();
	  shape.indices = (GLUSuint*)malloc(shape.numberIndices * sizeof(GLUSuint));
	  memcpy(shape.indices, &edgeIndices[0], sizeof(int) * shape.numberIndices);
  } else {
	  shape.numberIndices = triangles.size();
	  shape.indices = (GLUSuint*)malloc(shape.numberIndices * sizeof(GLUSuint));
	  memcpy(shape.indices, &triangles[0], sizeof(int) * shape.numberIndices);
  }
 
}

void ClothMeshData::updateNormalsGPU()
{
  if (pMesh == nullptr || this->vertNormals.size() == 0)
    return;
  GLUSshape& shape = pMesh->m_glusShape;
  shape.normals = (GLUSfloat*)malloc(3 * shape.numberVertices * sizeof(GLUSfloat));
  memcpy(shape.normals, &vertNormals[0], 3*sizeof(float) * shape.numberVertices);
  // copy current recalculated normals to appropriate VBO on GPU

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool isin(int x, int y, int el) {
	return (x >= 0) && (y >= 0) && (x < el) && (y < el);
}

void SimStep(ClothMeshData* pMesh, float delta_t)
{ 
	pMesh->global_time += delta_t;
	float g = ClothMeshData::g;
	float wind = ClothMeshData::wind;
	delta_t *= ClothMeshData::time_scale;
  // get in and out pointers
  //
  float4* inVertPos  = pMesh->pinPong ? &pMesh->vertPos1[0] : &pMesh->vertPos0[0];
  float4* inVertVel  = pMesh->pinPong ? &pMesh->vertVel1[0] : &pMesh->vertVel0[0];

  float4* outVertPos = pMesh->pinPong ? &pMesh->vertPos0[0] : &pMesh->vertPos1[0];
  float4* outVertVel = pMesh->pinPong ? &pMesh->vertVel0[0] : &pMesh->vertVel1[0];

  // accumulate forces first
  //
  for (size_t i = 0; i < pMesh->vertForces.size(); i++) // clear all forces
    pMesh->vertForces[i] = float4(0, 0, 0, 0);

  for (int connectId = 0; connectId < pMesh->connectionNumber(); connectId++)
  {
	  int a = pMesh->edgeIndices[2 * connectId];
	  int b = pMesh->edgeIndices[2 * connectId + 1];
	  float4 n = (inVertPos[b] - inVertPos[a]);
	  float dist = length(inVertPos[b] - inVertPos[a]);
	  float len = pMesh->edgeInitialLen[connectId];
	  if (dist < 1e-8)
		  continue;
	  n = normalize(n);
	  pMesh->vertForces[a] += n * pMesh->edgeHardness[connectId] * (dist - len);
	  pMesh->vertForces[b] -= n * pMesh->edgeHardness[connectId] * (dist - len);
  }
 

  // update positions and velocity
  //
  float el = ClothMeshData::elements;
  float ma = 0;
  for (int i = 0; i < pMesh->vertexNumber(); ++i) 
	  if (!pMesh->isFixed[i]) {
		  float wind2 = wind /**exp(-inVertPos[i].y)*/*exp(-pMesh->global_time/100)*cos(pMesh->global_time/100)*(1-abs(inVertPos[i].x));

		  float4 a = pMesh->vertForces[i] - float4(0, g, wind2, 0);
	      outVertPos[i] = inVertPos[i] + inVertVel[i] * delta_t + a * delta_t*delta_t / 2;
		  outVertVel[i] = ClothMeshData::loss*(inVertVel[i] + a * delta_t);
      }
  float planePos = 2.5 / 3;

  for (int i = 0; i < pMesh->vertexNumber(); ++i) 
	  if (outVertPos[i].y < -planePos){
		  //outVertPos[i].y = -planePos;
	  }

  int elem = ClothMeshData::elements;

  static const int dirs = 4;
  static int dx1[dirs] = { +1, 0, -1, 0 };
  static int dy1[dirs] = { 0, +1, 0, -1 };
  static int dx2[dirs] = { 0, -1, 0, +1 };
  static int dy2[dirs] = { +1, 0, -1, 0 };

  for (int i = 0; i < pMesh->vertexNumber(); ++i) {
	  int cnt = 0;
	  float4 sum = float4(0,0,0,0);

	  int x = i / ClothMeshData::elements;
	  int y = i % ClothMeshData::elements;

	  for (int j = 0; j < dirs; ++j) {
		  if (isin(x+dx1[j], y+dy1[j], elem) && isin(x + dx2[j], y + dy2[j], elem)) {
			  cnt++;
			  sum += cross(inVertPos[(x + dx1[j])*elem + y + dy1[j]] - inVertPos[x*elem + y],
				           inVertPos[(x + dx2[j])*elem + y + dy2[j]] - inVertPos[x*elem + y]);
		  }
	  }

	  sum /= cnt;
	  pMesh->vertNormals[i] = float3(sum.x, sum.y, sum.z);
  }
  pMesh->pinPong = !pMesh->pinPong; // swap pointers for next sim step
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void RecalculateNormals(ClothMeshData* pMesh)
{
}

