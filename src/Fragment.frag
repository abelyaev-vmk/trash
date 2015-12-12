#version 330

uniform sampler2D floor_texture;
uniform vec3 light_pos;
in vec2 fragmentTexCoord;
out vec4 fragColor;
in vec3 fragmentWorldPos;


void main(void)
{	
 // gl_Position *= 2;
  fragColor = texture(floor_texture, fragmentTexCoord*4);
  fragColor *= max(dot(normalize(light_pos-fragmentWorldPos), vec3(0,1,0)), 0.0f);
}





