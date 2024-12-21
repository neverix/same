#include <stdlib.h>
#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include <string.h>
#include <stdio.h>
#include "math.h"

#define maxi(a, b) ((a) > (b) ? (a) : (b))
#define imin(a, b) ((a) < (b) ? (a) : (b))
// #define BOARD_WIDTH 256
// #define BOARD_HEIGHT 256
#define BOARD_WIDTH 512
#define BOARD_HEIGHT 512
#define WALL_HEIGHT 5.0
#define GRAD_NOISE_BLOCK 16
#define GRAD_NOISE_THRESHOLD -0.1
#define ACCELERATION 0.5f
#define MOVE_SPEED 0.5f
#define BOUNCE_SPEED 0.5f

#define CAM_OFFSET 45.0
#define CAM_INITIAL_ANGLE (PI / 4)
#define MOUSE_SPEED 0.01f
#define MOUSE_X_SPEED 0.7f
#define MOUSE_Y_SPEED 0.3f

#define WORLD_RADIUS (BOARD_WIDTH / 2)
// #define WORLD_RADIUS 100.0
// #define WORLD_RADIUS 30.0
#define CYLINDER_RES 128
#define CYLINDER_HEIGHT (WALL_HEIGHT * 2)
#define CYLINDER_RADIUS (WORLD_RADIUS)
#define GROUND_EXTENTS 1000.0f

#define DESTRUCTION_RADIUS 5.0f
#define REGENERATE_MAX 10
#define REGEN_EVERY_S 2.0f

typedef struct Client {
    float width;
    float height;

    Camera3D camera;
    Vector2 camera_angle;
    bool camera_moved;

    Shader shader;

    Model outer_walls;
    Model walls;
} Client;

typedef enum Cell {
    EMPTY,
    WALL,
    RESOURCE,
} Cell;

typedef struct Entity {
    Vector2 position;
    Vector2 velocity;
    Vector2 target_velocity;
    float power;
    float atk;
    float def;
    float hp;
} Entity;

typedef struct {
    int n_entities;
    Entity *entities;

    Cell board[BOARD_HEIGHT][BOARD_WIDTH];
    bool was_destroyed[BOARD_HEIGHT][BOARD_WIDTH];
    float noise_buf[BOARD_HEIGHT][BOARD_WIDTH];
    int n_walls;
    int generation_seed;

    float elapsed_time;
    float last_regenned;
    int frame;
} Game;

void gen_board(Game *game);
void grad_noise(float *result_buf);

void draw(Client *client, Game *game);
void handle_controls(Client *client, Game *game);
bool process_collisions(Game *game, float dt);
bool regen_board(Game *game);

Mesh gen_cylinder();
void generate_wall(Client *client, Game *game);

inline float ent_radius(Entity *e) { return 0.5f; }

int main(void) {
    Game *game = calloc(1, sizeof(Game));
    game->n_entities = 50;
    game->entities = calloc(game->n_entities, sizeof(Entity));
    game->generation_seed = 15;

    gen_board(game);
    bool has_spawned[BOARD_HEIGHT][BOARD_WIDTH] = {0};
    for (int y = 0; y < BOARD_HEIGHT; y++) {
        for (int x = 0; x < BOARD_WIDTH; x++) {
            if (game->board[y][x] == WALL) {
                game->n_walls++;
                has_spawned[y][x] = true;
            }
        }
    }
    has_spawned[BOARD_HEIGHT / 2][BOARD_WIDTH / 2] = true;
    for (int i = 0; i < game->n_entities; i++) {
        for (int j = 0; j < 128; j++) {
            int y = rand() % BOARD_HEIGHT, x = rand() % BOARD_WIDTH;
            Vector2 pos = {x - BOARD_WIDTH / 2, y - BOARD_HEIGHT / 2};
            float radius = ent_radius(&game->entities[i]);
            if (Vector2Length(pos) < WORLD_RADIUS - radius && !has_spawned[y][x]) {
                game->entities[i] = (Entity) {0};
                game->entities[i].position = pos;
                has_spawned[y][x] = true;
                break;
            }
        }
    }

    Client *client = calloc(1, sizeof(Client));
    client->width = 800;
    client->height = 600;
    client->camera = (Camera3D){(Vector3){}, (Vector3){},
                                (Vector3){0.0f, 1.0f, 0.0f}, 45.0f, 0};
    client->camera_angle = (Vector2){CAM_INITIAL_ANGLE, 0.0f};

    InitWindow(client->width, client->height, "same");

    Shader shader = LoadShader("vert.vs", "frag.fs");

    client->shader = shader;

    Mesh cylinder = gen_cylinder();
    client->outer_walls = LoadModelFromMesh(cylinder);
    client->outer_walls.materials[0].shader = client->shader;

    generate_wall(client, game);

    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        float dt = GetFrameTime();
        draw(client, game);
        handle_controls(client, game);
        if (game->frame % 120 == 0) {
            for (int i = 1; i < game->n_entities; i++) {
                Vector2 target = (Vector2) {rand() / (float)RAND_MAX,
                                            rand() / (float)RAND_MAX};
                target = Vector2Scale(Vector2Normalize(target), MOVE_SPEED);
                game->entities[i].target_velocity = target;
            }
        }
        if(process_collisions(game, dt) == true) {
            generate_wall(client, game);
        }
        if (game->elapsed_time - game->last_regenned > REGEN_EVERY_S) {
            if (regen_board(game)) {
                generate_wall(client, game);
            }
            game->last_regenned = game->elapsed_time;
        }
        game->frame++;
    }

    CloseWindow();
    free(client);
    return 0;
}

bool process_collisions(Game *game, float dt) {
    bool wall_changed = false;
    for (int i = 0; i < game->n_entities; i++) {
        game->entities[i].velocity = Vector2Add(
            game->entities[i].velocity,
            Vector2ClampValue(Vector2Subtract(game->entities[i].target_velocity,
                                              game->entities[i].velocity),
                              0.0, ACCELERATION * dt));
        Vector2 og_velocity = game->entities[i].velocity;
        Vector2 original_position = game->entities[i].position;
        Vector2 projected_position = Vector2Add(original_position, og_velocity);
        bool collided = false;
        bool with_wall = false;
        Vector2 collision_direction;

        float radius = ent_radius(&game->entities[i]);
        if (Vector2Length(projected_position) > WORLD_RADIUS - radius) {
            collided = true;
            collision_direction = Vector2Negate(Vector2Normalize(projected_position));
        } else {
            int consider_dirs = maxi(10, (int)radius * PI * 3);
            Vector2 directions = {0, 0};
            int flagged = 0;
            for (int dir = 0; dir < consider_dirs; dir++) {
                float angle = dir / (float)consider_dirs * 2 * PI;
                Vector2 extrapolation = {cosf(angle), sinf(angle)};
                extrapolation =
                    Vector2Scale(Vector2Normalize(extrapolation), radius);
                Vector2 extrapolated_position =
                    Vector2Add(projected_position, extrapolation);
                if (game->board[(int)(extrapolated_position.y +
                                      BOARD_HEIGHT / 2)][(
                        int)(extrapolated_position.x + BOARD_WIDTH / 2)] ==
                    WALL) {
                    flagged++;
                    directions = Vector2Add(directions, extrapolation);
                }
            }
            if (flagged > 0) {
                collided = true;
                with_wall = true;
                collision_direction =
                    Vector2Scale(Vector2Normalize(directions), -1.0f);
            }
        }
        if (collided) {
            game->entities[i].velocity = Vector2Subtract(
                og_velocity, Vector2Scale(collision_direction,
                                          Vector2DotProduct(collision_direction,
                                                            og_velocity)));
            game->entities[i].velocity =
                Vector2Add(game->entities[i].velocity,
                           Vector2Scale(collision_direction, BOUNCE_SPEED));
            Vector2 norm_velocity = Vector2Normalize(og_velocity);
            game->entities[i].target_velocity = Vector2Subtract(
                game->entities[i].target_velocity,
                Vector2Scale(
                    norm_velocity,
                    Vector2DotProduct(norm_velocity,
                                      game->entities[i].target_velocity)));
            game->entities[i].position = original_position;
        } else {
            game->entities[i].position = projected_position;
        }

        if (with_wall) {
            Vector2 epicenter = game->entities[i].position;
            for (int ox = -ceil(DESTRUCTION_RADIUS); ox <= ceil(DESTRUCTION_RADIUS); ox++) {
                for (int oy = -ceil(DESTRUCTION_RADIUS); oy <= ceil(DESTRUCTION_RADIUS); oy++) {
                    Vector2 offset = {ox, oy};
                    if (Vector2Length(offset) > DESTRUCTION_RADIUS) {
                        continue;
                    }
                    Vector2 target = Vector2Add(epicenter, offset);
                    if (target.x < -BOARD_WIDTH / 2 || target.x >= BOARD_WIDTH / 2 ||
                        target.y < -BOARD_HEIGHT / 2 || target.y >= BOARD_HEIGHT / 2) {
                        continue;
                    }
                    if (game->board[(int)(target.y + BOARD_HEIGHT / 2)][(int)(target.x + BOARD_WIDTH / 2)] == WALL) {
                        int a = (int)(target.y + BOARD_HEIGHT / 2), b = (int)(target.x + BOARD_WIDTH / 2);
                        game->board[a][b] = EMPTY;
                        game->was_destroyed[a][b] = true;
                        wall_changed = true;
                    }
                }
            }
        }
    }
    game->elapsed_time += dt;
    return wall_changed;
}

void handle_controls(Client *client, Game *game) {
    Vector2 camera_motion = GetMouseDelta();
    if (Vector2Length(camera_motion) > 0.0) {
        if (client->camera_moved)
            camera_motion = (Vector2){camera_motion.y * MOUSE_Y_SPEED,
                                      camera_motion.x * MOUSE_X_SPEED};
        else
            camera_motion = (Vector2){0.0, 0.0};
        Vector2 new_camera = Vector2Add(
            client->camera_angle, Vector2Scale(camera_motion, MOUSE_SPEED));
        new_camera.x = Clamp(new_camera.x, 0.0, PI * 0.4);
        client->camera_angle = new_camera;
        client->camera_moved = true;
    }

    game->entities[0].target_velocity = (Vector2){0, 0};

    Vector2 camera_forward = Vector2Normalize(Vector2Subtract(
        (Vector2){client->camera.target.x, client->camera.target.z},
        (Vector2){client->camera.position.x, client->camera.position.z}));
    camera_forward = Vector2Scale(camera_forward, MOVE_SPEED);
    Vector2 camera_right = Vector2Rotate(camera_forward, PI / 2);
    if (IsKeyDown(KEY_K)) {
        game->entities[0].target_velocity =
            Vector2Add(game->entities[0].target_velocity, camera_forward);
    }
    if (IsKeyDown(KEY_J)) {
        game->entities[0].target_velocity =
            Vector2Subtract(game->entities[0].target_velocity, camera_forward);
    }
    if (IsKeyDown(KEY_L)) {
        game->entities[0].target_velocity =
            Vector2Add(game->entities[0].target_velocity, camera_right);
    }
    if (IsKeyDown(KEY_H)) {
        game->entities[0].target_velocity =
            Vector2Subtract(game->entities[0].target_velocity, camera_right);
    }
    game->entities[0].target_velocity =
        Vector2ClampValue(game->entities[0].target_velocity, 0.0, MOVE_SPEED);
}

void draw(Client *client, Game *game) {
    BeginDrawing();
    ClearBackground(DARKGRAY);

    BeginMode3D(client->camera);

    // DrawGrid(maxi(BOARD_WIDTH, BOARD_HEIGHT), 1.0f);


    BeginShaderMode(client->shader);

    DrawModel(client->outer_walls, (Vector3){0.0f, 0.0f, 0.0f}, 1.0f,
              (Color){255, 255, 0, 255});

    DrawPlane((Vector3){0.0f, 0.0f, 0.0f}, (Vector2){GROUND_EXTENTS, GROUND_EXTENTS}, (Color){0, 255, 0, 255});

    Material mat = LoadMaterialDefault();
    mat.maps[MATERIAL_MAP_DIFFUSE].color = (Color){0, 255, 255, 0};
    DrawModel(client->walls, (Vector3){-BOARD_WIDTH / 2 + 1, 0.0f, -BOARD_HEIGHT / 2 + 1}, 1.0f, WHITE);

    for (int i = 0; i < game->n_entities; i++) {
        float radius = ent_radius(&game->entities[i]);
        DrawSphere((Vector3){game->entities[i].position.x, radius,
                             game->entities[i].position.y},
                   radius, RED);
    }


    client->camera.target = (Vector3){game->entities[0].position.x, 0.0f,
                                      game->entities[0].position.y};
    Vector3 camera_position = {0, 0, CAM_OFFSET};
    Vector2 new_yz = Vector2Rotate((Vector2){camera_position.y, camera_position.z}, -client->camera_angle.x);
    camera_position.y = new_yz.x;
    camera_position.z = new_yz.y;
    Vector2 new_xz = Vector2Rotate((Vector2){camera_position.x, camera_position.z}, -client->camera_angle.y);
    camera_position.x = new_xz.x;
    camera_position.z = new_xz.y;
    client->camera.position = Vector3Add(
        (Vector3){game->entities[0].position.x, 0.0,
                  game->entities[0].position.y},
        camera_position);

    EndShaderMode();

    EndMode3D();

    DrawFPS(10, 10);
    char textbuf[256];
    snprintf(textbuf, 256,
             "x: %.1f, y: %.1f, vx: %.1f, vy: %.1f, tvx: %.1f, tvy: %.1f",
             game->entities[0].position.x, game->entities[0].position.y,
             game->entities[0].velocity.x, game->entities[0].velocity.y,
             game->entities[0].target_velocity.x,
             game->entities[0].target_velocity.y);
    DrawText(textbuf, 100, 10, 20, RED);

    EndDrawing();
}

void grad_noise(float *result_buf) {
    float grad_noise_base[BOARD_HEIGHT / GRAD_NOISE_BLOCK]
                         [BOARD_WIDTH / GRAD_NOISE_BLOCK][2];
    for (int y = 0; y < BOARD_HEIGHT / GRAD_NOISE_BLOCK; y++) {
        for (int x = 0; x < BOARD_WIDTH / GRAD_NOISE_BLOCK; x++) {
            float u = (float)rand() / RAND_MAX - 0.5,
                  v = (float)rand() / RAND_MAX - 0.5;
            float len = sqrt(u * u + v * v);
            grad_noise_base[y][x][0] = u / len;
            grad_noise_base[y][x][1] = v / len;
        }
    }
    float grad_noise[BOARD_HEIGHT][BOARD_WIDTH];
    for (int y = 0; y < BOARD_HEIGHT; y += GRAD_NOISE_BLOCK) {
        for (int x = 0; x < BOARD_WIDTH; x += GRAD_NOISE_BLOCK) {
            float cu =
                grad_noise_base[y / GRAD_NOISE_BLOCK][x / GRAD_NOISE_BLOCK][0];
            float cv =
                grad_noise_base[y / GRAD_NOISE_BLOCK][x / GRAD_NOISE_BLOCK][1];
            for (int i = 0; i < GRAD_NOISE_BLOCK; i++) {
                for (int j = 0; j < GRAD_NOISE_BLOCK; j++) {
                    int ry = y + i;
                    int rx = x + j;
                    float dy = i / (float)GRAD_NOISE_BLOCK * 2 - 1;
                    float dx = j / (float)GRAD_NOISE_BLOCK * 2 - 1;
                    float current_v = dy * cu + dx * cv;
                    grad_noise[ry][rx] = current_v;
                }
            }
        }
    }

    for (int y = 0; y < BOARD_HEIGHT; y++) {
        for (int x = 0; x < BOARD_WIDTH; x++) {
            float accum = 0;
            int total = 0;
            for (int yo = -GRAD_NOISE_BLOCK / 2; yo < GRAD_NOISE_BLOCK / 2;
                 yo++) {
                for (int xo = -GRAD_NOISE_BLOCK / 2; xo < GRAD_NOISE_BLOCK / 2;
                     xo++) {
                    int ry = y + yo;
                    int rx = x + xo;
                    if (ry < 0 || ry >= BOARD_HEIGHT || rx < 0 ||
                        rx >= BOARD_WIDTH) {
                        continue;
                    }
                    accum += grad_noise[ry][rx];
                    total++;
                }
            }
            float current_v = accum / total;
            result_buf[y * BOARD_WIDTH + x] = current_v;
        }
    }
}

inline void gen_board(Game *game) {
    srand(game->generation_seed);
    grad_noise((float*)game->noise_buf);
    for (int y = 0; y < BOARD_HEIGHT; y++) {
        for (int x = 0; x < BOARD_WIDTH; x++) {
            Vector2 pos = {x - BOARD_WIDTH / 2, y - BOARD_HEIGHT / 2};
            pos = Vector2Add(pos, Vector2Normalize(pos));
            if (Vector2Length(pos) > WORLD_RADIUS) {
                game->board[y][x] = EMPTY;
                continue;
            }
            if (game->noise_buf[y][x] > GRAD_NOISE_THRESHOLD) {
                game->board[y][x] = EMPTY;
            } else {
                game->board[y][x] = WALL;
            }
        }
    }
}

int cmpfunc(const void *a, const void *b) {
    float x = *(float*)a;
    float y = *(float*)b;
    if (x > y) {
        return 1;
    } else if (x < y) {
        return -1;
    } else {
        return 0;
    }
}

inline bool regen_board(Game *game) {
    int current_walls = 0;
    for (int y = 0; y < BOARD_HEIGHT; y++) {
        for (int x = 0; x < BOARD_WIDTH; x++) {
            if (game->board[y][x] == WALL) {
                current_walls++;
            }
        }
    }
    int to_generate = imin(game->n_walls - current_walls, REGENERATE_MAX * game->n_entities);
    // printf("to_generate: %d\n", to_generate);
    if (to_generate == 0) {
        return false;
    }

    float eligible_noises[BOARD_HEIGHT * BOARD_WIDTH] = {__FLT_MAX__};
    int n_eligible = 0;

    bool players_present[BOARD_HEIGHT][BOARD_WIDTH] = {0};
    for (int i = 0; i < game->n_entities; i++) {
        Vector2 pos = game->entities[i].position;
        int y = pos.y + BOARD_HEIGHT / 2, x = pos.x + BOARD_WIDTH / 2;
        for (int oy = -2; oy <= 2; oy++) {
            for (int ox = -2; ox <= 2; ox++) {
                int ry = y + oy, rx = x + ox;
                if (ry < 0 || ry >= BOARD_HEIGHT || rx < 0 || rx >= BOARD_WIDTH) {
                    continue;
                }
                players_present[ry][rx] = true;
            }
        }
    }

    #define is_eligible(y, x) (game->was_destroyed[y][x] && game->board[y][x] == EMPTY \
        && Vector2Length((Vector2){(x) - BOARD_WIDTH / 2, (y) - BOARD_HEIGHT / 2}) < WORLD_RADIUS \
        && !players_present[y][x])

    for (int y = 0; y < BOARD_HEIGHT; y++) {
        for (int x = 0; x < BOARD_WIDTH; x++) {
            if (is_eligible(y, x)) {
                eligible_noises[n_eligible++] = game->noise_buf[y][x];
            }
        }
    }
    if (n_eligible < to_generate) {
        to_generate = n_eligible;
    }
    // printf("eligible: %d\n", to_generate);
    qsort(eligible_noises, BOARD_HEIGHT * BOARD_WIDTH, sizeof(float), cmpfunc);
    float regen_cutoff = eligible_noises[to_generate - 1];
    bool changed = false;
    for (int y = 0; y < BOARD_HEIGHT; y++) {
        for (int x = 0; x < BOARD_WIDTH; x++) {
            if (is_eligible(y, x) && game->noise_buf[y][x] <= regen_cutoff) {
                game->board[y][x] = WALL;
                game->was_destroyed[y][x] = false;
                changed = true;
            }
        }
    }
    return changed;
}

inline Mesh gen_cylinder() {
    Mesh cylinder = {0};
    cylinder.triangleCount = CYLINDER_RES * 2;
    cylinder.vertexCount = cylinder.triangleCount * 3;
    cylinder.vertices = (float *)MemAlloc(
        cylinder.vertexCount * 3 *
        sizeof(float));  // 3 vertices, 3 coordinates each (x, y, z)
    cylinder.normals = (float *)MemAlloc(
        cylinder.vertexCount * 3 *
        sizeof(float));  // 3 vertices, 3 coordinates each (x, y, z)
    cylinder.texcoords = (float *)MemAlloc(
        cylinder.vertexCount * 2 *
        sizeof(float));  // 3 vertices, 2 coordinates each (x, y)

    for (int quad = 0; quad < CYLINDER_RES; quad++) {
        Vector3 bottom_left = {
            cosf(quad / (float)CYLINDER_RES * 2 * PI) * CYLINDER_RADIUS, 0.0,
            sinf(quad / (float)CYLINDER_RES * 2 * PI) * CYLINDER_RADIUS};
        Vector3 bottom_right = {
            cosf((quad - 1) / (float)CYLINDER_RES * 2 * PI) * CYLINDER_RADIUS,
            0.0,
            sinf((quad - 1) / (float)CYLINDER_RES * 2 * PI) * CYLINDER_RADIUS};
        Vector3 top_left = {bottom_left.x, CYLINDER_HEIGHT, bottom_left.z};
        Vector3 top_right = {bottom_right.x, CYLINDER_HEIGHT, bottom_right.z};
        Vector3 normal = Vector3Normalize(
            Vector3CrossProduct(Vector3Subtract(bottom_left, top_left),
                                Vector3Subtract(bottom_right, bottom_left)));
        int base = quad * 18;
        for (int i = 0; i < 6; i++) {
            cylinder.normals[base + i * 3] = normal.x;
            cylinder.normals[base + i * 3 + 1] = normal.y;
            cylinder.normals[base + i * 3 + 2] = normal.z;
        }

        cylinder.vertices[base] = bottom_left.x;
        cylinder.vertices[base + 1] = bottom_left.y;
        cylinder.vertices[base + 2] = bottom_left.z;
        cylinder.vertices[base + 3] = top_left.x;
        cylinder.vertices[base + 4] = top_left.y;
        cylinder.vertices[base + 5] = top_left.z;
        cylinder.vertices[base + 6] = top_right.x;
        cylinder.vertices[base + 7] = top_right.y;
        cylinder.vertices[base + 8] = top_right.z;
        cylinder.vertices[base + 9] = top_right.x;
        cylinder.vertices[base + 10] = top_right.y;
        cylinder.vertices[base + 11] = top_right.z;
        cylinder.vertices[base + 12] = bottom_right.x;
        cylinder.vertices[base + 13] = bottom_right.y;
        cylinder.vertices[base + 14] = bottom_right.z;
        cylinder.vertices[base + 15] = bottom_left.x;
        cylinder.vertices[base + 16] = bottom_left.y;
        cylinder.vertices[base + 17] = bottom_left.z;
    }
    UploadMesh(&cylinder, false);
    return cylinder;
}

inline void generate_wall(Client *client, Game *game) {
    if (client->walls.meshes != NULL) {
        UnloadModel(client->walls);
    }
    unsigned char wall_data[BOARD_HEIGHT][BOARD_WIDTH][3] = {0};
    Image walls = (Image){.width = BOARD_WIDTH,
                          .height = BOARD_HEIGHT,
                          .format = PIXELFORMAT_UNCOMPRESSED_R8G8B8,
                          .data = (void *)wall_data};
    for (int y = 0; y < BOARD_HEIGHT; y++) {
        for (int x = 0; x < BOARD_WIDTH; x++) {
            unsigned char col = game->board[y][x] == WALL ? 255 : 0;
            wall_data[y][x][0] = col;
            wall_data[y][x][1] = col;
            wall_data[y][x][2] = col;
        }
    }
    Mesh walls_mesh = GenMeshCubicmap(walls, (Vector3){1.0, WALL_HEIGHT, 1.0});
    client->walls = LoadModelFromMesh(walls_mesh);
    client->walls.materials[0].shader = client->shader;
    client->walls.materials[0].maps[MATERIAL_MAP_DIFFUSE].color =
        (Color){0, 255, 255, 255};
}