#version 330

uniform sampler2D floor_texture;

in vec3 fragmentWorldPos;
in vec3 fragmentNormal;
in vec2 fragmentTexCoord;

in vec3 fragmentSunView;
in vec3 fragmentNormalView;

uniform int render_ind;
uniform vec3 light_pos;

out vec4 fragColor;

void main(void)
{	

  vec3 light = normalize(light_pos - fragmentWorldPos);
  if (render_ind == 2)
      fragColor = vec4(fragmentNormal, 1);
  else if (render_ind == 1)
      fragColor = vec4(1, 0, 0, 1);
  else {
      fragColor = texture(floor_texture, fragmentTexCoord);
	  fragColor *= max(dot(light, fragmentNormal), 0.0f);
  }
}



