
#version 330

// Input vertex attributes (from vertex shader)
in vec2 fragTexCoord;
in vec4 fragColor;
in vec3 fragPosition;
in vec3 fragPositionWorld;
in vec3 fragNormal;
in vec3 rayDir;

uniform mat4 mvp;

uniform vec4 colDiffuse;

// Output fragment color
out vec4 finalColor;

// float max_terrain_height = 1.0;
float max_terrain_height = 2.0;

#define fbm_steps 4
// #define raymarch_steps 2
#define raymarch_steps 8
#define wall_height 5.0
#define outer_wall_height 10.0

vec4 gray(float color) {
    return vec4(color, color, color, 1.0);
}

float rand(vec2 co) {
    return fract(sin(dot(co, vec2(839.27, -312.39))) * 21890.29);
}

float value_noise(vec2 co) {
    vec2 base = floor(co);
    float ul = rand(base);
    float ur = rand(base + vec2(1.0, 0.0));
    float ll = rand(base + vec2(0.0, 1.0));
    float lr = rand(base + vec2(1.0, 1.0));
    float u = mix(ul, ur, fract(co.x));
    float l = mix(ll, lr, fract(co.x));
    return mix(u, l, fract(co.y));
}

float rand_3d(vec3 co) {
    return fract(sin(dot(co, vec3(839.27, -312.39, 895.12))) * 21890.29);
}

float value_noise_3d(vec3 co) {
    vec3 base = floor(co);
    float ulb = rand_3d(base);
    float ulf = rand_3d(base + vec3(1.0, 0.0, 0.0));
    float urb = rand_3d(base + vec3(0.0, 1.0, 0.0));
    float urf = rand_3d(base + vec3(1.0, 1.0, 0.0));
    float llb = rand_3d(base + vec3(0.0, 0.0, 1.0));
    float llf = rand_3d(base + vec3(1.0, 0.0, 1.0));
    float lrb = rand_3d(base + vec3(0.0, 1.0, 1.0));
    float lrf = rand_3d(base + vec3(1.0, 1.0, 1.0));
    float ul = mix(ulb, ulf, fract(co.x));
    float ur = mix(urb, urf, fract(co.x));
    float ll = mix(llb, llf, fract(co.x));
    float lr = mix(lrb, lrf, fract(co.x));
    float u = mix(ul, ur, fract(co.y));
    float l = mix(ll, lr, fract(co.y));
    return mix(u, l, fract(co.z));
}

float fbm(vec2 st) {
    vec2 coord = st;
    float value = 0.0;
    float amp_mul = 0.5;
    float amplitude = 1.0 - amp_mul;
    float frequency = 1.0;
    vec2 offset = vec2(0.312, -0.732);
    float angle = 0.13;
    // mat2 rotation = mat2(cos(angle), sin(angle), -sin(angle), cos(angle));
    for (int i = 0; i < fbm_steps; i++) {
        value += amplitude * value_noise(coord * frequency);
        coord += offset;
        // coord = rotation * coord;
        frequency *= 2.0;
        amplitude *= amp_mul;
    }
    return value;
}

float fbm_3d(vec3 st) {
    vec3 coord = st;
    float value = 0.0;
    float amp_mul = 0.5;
    float amplitude = 1.0 - amp_mul;
    float frequency = 1.0;
    // vec3 offset = vec3(0.312, -0.732, 0.183);
    vec3 offset = vec3(-0.732, 0.312, 0.183);
    float angle = 0.13;
    // mat2 rotation = mat2(cos(angle), sin(angle), -sin(angle), cos(angle));
    for (int i = 0; i < fbm_steps; i++) {
        value += amplitude * value_noise_3d(coord * frequency);
        coord += offset;
        // coord = rotation * coord;
        // coord.xz = rotation * coord.xz;
        // coord.xy = rotation * coord.xy;
        frequency *= 2.0;
        amplitude *= amp_mul;
    }
    return value;
}

float terrain_height(vec2 st) {
    return max_terrain_height * fbm(st * 0.2);
}

vec3 darkGreen = vec3(0.1, 0.2, 0.1) / 0.5;
vec3 slightlyLighterGreen = vec3(0.2, 0.3, 0.2) / 0.5;
vec3 darkBrown = vec3(0.2, 0.1, 0.0) / 0.6;
vec3 brown = vec3(0.35, 0.25, 0.1) / 0.6;
vec3 darkestBrown = darkBrown * (darkBrown / brown);

void main() {
    if (fragColor.r == 0.0 && fragColor.g == 0.0 && fragColor.b == 0.0) {
        // finalColor = gray(fragPositionWorld.y / 70.0 + 0.1);
        float height = fragPositionWorld.y / wall_height;
        finalColor = vec4(mix(darkBrown, brown, sqrt(height)), 1.0);
    } else if (fragColor.r == 0.0 && fragColor.g == 1.0 && fragColor.b == 0.0) {
        // vec3 rayDir = normalize(inverse(mat3(mvp)) * vec3(0.0, 0.0, 1.0));
        float rayLength = max_terrain_height / (-rayDir.y);
        vec3 rayOrigin = fragPositionWorld - rayDir * rayLength - vec3(0.0, max_terrain_height, 0.0);

        for (int i = 0; i < raymarch_steps; i++) {
            float noiseValue = terrain_height(rayOrigin.xz);
            if (rayOrigin.y + max_terrain_height < noiseValue) {
                break;
            }
            rayOrigin += rayDir * rayLength / raymarch_steps;
        }

        float noiseValue = terrain_height(rayOrigin.xz) / max_terrain_height;
        vec3 color = mix(darkestBrown, darkBrown, noiseValue);

        finalColor = vec4(color, 1.0);
    } else if (colDiffuse.r == 1.0 && colDiffuse.g == 1.0 && colDiffuse.b == 0.0) {
        // vec3 index = fragPositionWorld / 10.0;
        // float x = fbm_3d(index);
        // float y = fbm_3d(index + vec3(0.0, 0.0, 10.0));
        // float z = fbm_3d(index + vec3(0.0, 0.0, 2-.0));
        // vec3 new_xyz = index + vec3(x, y, z) * 10;
        // float val = fbm_3d(new_xyz);
        
        // // finalColor = gray(val);
        // // vec3 darkPurple = vec3(0.1, 0.0, 0.1);
        // // vec3 lightPurple = vec3(0.5, 0.0, 0.7);
        // vec3 darkPurple = darkestBrown;
        // // vec3 lightPurple = brown;
        // vec3 lightPurple = vec3(0.5, 0.0, 0.7);
        // finalColor = vec4(mix(darkPurple, lightPurple, val), 1.0);

        float height = fragPositionWorld.y / outer_wall_height;
        finalColor = vec4(mix(darkBrown, brown, height), 1.0);
    } else {
        finalColor = fragColor;
    }
}
