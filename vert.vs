
#version 330

// Input vertex attributes
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;
in vec4 vertexColor;

// Input uniform values
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;

// Output vertex attributes (to fragment shader)
out vec3 fragPosition;
out vec3 fragPositionWorld;
out vec2 fragTexCoord;
out vec4 fragColor;
out vec3 fragNormal;
out vec3 rayDir;

// NOTE: Add here your custom variables

void main() {
    // Send vertex attributes to fragment shader
    fragPositionWorld = vertexPosition;
    fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));
    fragTexCoord = vertexTexCoord;
    fragColor = vertexColor;
    // fragNormal = normalize(vec3(matNormal * vec4(vertexNormal, 1.0)));
    fragNormal = normalize(vertexNormal);

    // Calculate final vertex position
    gl_Position = mvp * vec4(vertexPosition, 1.0);
    
    rayDir = normalize(inverse(mat3(mvp)) * vec3(0.0, 0.0, 1.0));
}
