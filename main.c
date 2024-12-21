#include <stdlib.h>
#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include <string.h>
#include <stdio.h>
#include "math.h"

#define maxi(a, b) ((a) > (b) ? (a) : (b))
#define BOARD_WIDTH 256
#define BOARD_HEIGHT 256
#define WALL_HEIGHT 5.0
#define GRAD_NOISE_BLOCK 16
#define GRAD_NOISE_THRESHOLD -0.1
#define ACCELERATION 0.5f
#define MOVE_SPEED 0.3f
#define BOUNCE_SPEED 0.3f

#define CAM_OFFSET 45.0
#define CAM_INITIAL_ANGLE (PI / 4)
#define MOUSE_SPEED 0.01f
#define MOUSE_X_SPEED 0.7f
#define MOUSE_Y_SPEED 0.3f

// #define WORLD_RADIUS 100.0
#define WORLD_RADIUS 30.0
#define CYLINDER_RES 128
#define CYLINDER_HEIGHT (WALL_HEIGHT * 2)
#define CYLINDER_RADIUS (WORLD_RADIUS)
#define GROUND_EXTENTS 1000.0f

typedef struct Client {
    float width;
    float height;

    Camera3D camera;
    Vector2 camera_angle;
    bool camera_moved;

    Shader shader;

    Model outer_walls;
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
    Cell board[BOARD_HEIGHT][BOARD_WIDTH];
    int n_entities;
    Entity *entities;
} Game;

void gen_board(Game *game);
void grad_noise(float *result_buf);

void draw(Client *client, Game *game);
void handle_controls(Client *client, Game *game);
void process_collisions(Game *game, float dt);

Mesh gen_cylinder();

inline float ent_radius(Entity *e) { return 0.5f; }

int main(void) {
    Game *game = calloc(1, sizeof(Game));
    game->n_entities = 1;
    game->entities = calloc(game->n_entities, sizeof(Entity));
    game->entities[0].position = (Vector2){0, 0};

    gen_board(game);

    Client *client = calloc(1, sizeof(Client));
    client->width = 800;
    client->height = 600;
    client->camera = (Camera3D){(Vector3){}, (Vector3){},
                                (Vector3){0.0f, 1.0f, 0.0f}, 45.0f, 0};
    client->camera_angle = (Vector2){CAM_INITIAL_ANGLE, 0.0f};

    InitWindow(client->width, client->height, "same");

    client->shader = LoadShader("vert.vs", "frag.fs");

    Mesh cylinder = gen_cylinder();
    client->outer_walls = LoadModelFromMesh(cylinder);
    client->outer_walls.materials[0].shader = client->shader;

    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        float dt = GetFrameTime();
        draw(client, game);
        handle_controls(client, game);
        process_collisions(game, dt);
    }

    CloseWindow();
    free(client);
    return 0;
}

void process_collisions(Game *game, float dt) {
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
    }
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

    for (int i = 0; i < game->n_entities; i++) {
        game->entities[i].target_velocity = (Vector2){0, 0};
    }

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

    for (int y = 0; y < BOARD_HEIGHT; y++) {
        for (int x = 0; x < BOARD_WIDTH; x++) {
            if (game->board[y][x] == WALL) {
                DrawCube((Vector3){x - BOARD_WIDTH / 2 + 0.5f, WALL_HEIGHT / 2,
                                   y - BOARD_HEIGHT / 2 + 0.5f},
                         1.0f, WALL_HEIGHT, 1.0f, BLACK);
            }
        }
    }

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

void gen_board(Game *game) {
    float noise_buf[BOARD_HEIGHT][BOARD_WIDTH];
    grad_noise((float *)noise_buf);
    for (int y = 0; y < BOARD_HEIGHT; y++) {
        for (int x = 0; x < BOARD_WIDTH; x++) {
            Vector2 pos = {x - BOARD_WIDTH / 2, y - BOARD_HEIGHT / 2};
            pos = Vector2Add(pos, Vector2Normalize(pos));
            if (Vector2Length(pos) > WORLD_RADIUS) {
                game->board[y][x] = EMPTY;
                continue;
            }
            if (noise_buf[y][x] > GRAD_NOISE_THRESHOLD) {
                game->board[y][x] = EMPTY;
            } else {
                game->board[y][x] = WALL;
            }
        }
    }
}

Mesh gen_cylinder() {
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