#include <stdlib.h>
#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include <string.h>
#include <stdio.h>
#include "math.h"

#define maxi(a, b) ((a) > (b) ? (a) : (b))
#define imin(a, b) ((a) < (b) ? (a) : (b))
#define sign(x) ((x) > 0 ? 1 : -1)
// #define BOARD_WIDTH 256
// #define BOARD_HEIGHT 256
#define BOARD_WIDTH 512
#define BOARD_HEIGHT 512
// #define BOARD_WIDTH 1024
// #define BOARD_HEIGHT 1024
#define WALL_HEIGHT 5.0f
#define RENDER_CHUNK 32
#define GRAD_NOISE_BLOCK 16
// #define GRAD_NOISE_THRESHOLD -0.05
#define GRAD_NOISE_THRESHOLD -0.01
// #define ACCELERATION 0.5f
#define ACCELERATION 2.0f
#define MOVE_SPEED 0.5f
#define ATTRACTION_MULTIPLIER_VELOCITY 2.0f
#define ATTRACTION_MULTIPLIER_ACCELERATION 3.0f
#define BOUNCE_SPEED 0.2f

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

#define SKY_SIZE 1000.0f

#define MIN_INIT_SIZE 0.2f
#define MAX_INIT_SIZE 0.5f
#define MIN_RADIUS 1.5f
// #define MIN_RADIUS 2.0f

#define LOG_SIZE_DECAY_PER_S 0.03f
#define LOG_ATTRACTION_DECAY_PER_S 0.3f
/* #define LOG_ATTRACTION_DECAY_PER_S 1.2f */
#define PER_DESTROYED_SIZE 0.0004f
// #define LOSE_BOUNCE_BONUS 0.03f
#define LOSE_BOUNCE_BONUS 0.0f
#define BALL_TO_BALL_IMPULSE (MOVE_SPEED * 2)
#define ATTRACT_SIZE 0.25f

#define DEATH_SIZE 0.1f

#define FIGHT_PERIOD 0.2f

typedef enum Cell {
    EMPTY,
    WALL,
    RESOURCE,
} Cell;

typedef struct Client {
    float width;
    float height;

    Camera3D camera;
    Vector2 camera_angle;
    bool camera_moved;

    Shader shader;
    Mesh sky_wall;
    Model outer_walls;
    Model walls[BOARD_HEIGHT / RENDER_CHUNK][BOARD_WIDTH / RENDER_CHUNK];
    Cell past_board[BOARD_HEIGHT][BOARD_WIDTH];
} Client;

typedef struct Entity {
    Vector2 position;
    Vector2 velocity;
    Vector2 target_velocity;
    float size;
    float immune_since;
    bool attracted;
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

void destroy_entity(Game *game, Entity *e);

void draw(Client *client, Game *game);
void handle_controls(Client *client, Game *game);
bool process_collisions(Game *game, float dt);
bool regen_board(Game *game);
void copy_past_board(Client *client, Game *game);
void attract_to_nearest(Entity *e, Game *game);

Mesh gen_cylinder();
void generate_wall(Client *client, Game *game);

inline float ent_radius(Entity *e) {
    return sqrt(e->size / MIN_INIT_SIZE) * MIN_RADIUS;
}

inline float deceleration(Entity *e) {
    if (e-> size > MIN_INIT_SIZE)
        return e->size / MIN_INIT_SIZE;
    else
        return 1.0;
}

int main(void) {
    Game *game = calloc(1, sizeof(Game));
    game->n_entities = 128;
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
        game->entities[i].size =
            MIN_INIT_SIZE +
            (MAX_INIT_SIZE - MIN_INIT_SIZE) * rand() / (float)RAND_MAX;
    }

    Client *client = calloc(1, sizeof(Client));
    client->width = 800;
    client->height = 600;
    client->camera = (Camera3D){(Vector3){}, (Vector3){},
                                (Vector3){0.0f, 1.0f, 0.0f}, 45.0f, 0};
    client->camera_angle = (Vector2){CAM_INITIAL_ANGLE, 0.0f};

    InitWindow(client->width, client->height, "same");

    Shader shader = LoadShader("vert.glsl", "frag.glsl");

    client->shader = shader;

    Mesh cylinder = gen_cylinder();
    client->outer_walls = LoadModelFromMesh(cylinder);
    client->outer_walls.materials[0].shader = client->shader;

    Mesh plane = GenMeshPlane(SKY_SIZE, SKY_SIZE, 1, 1);
    client->sky_wall = plane;

    generate_wall(client, game);

    SetTargetFPS(60);
    SetTraceLogLevel(LOG_WARNING);

    while (!WindowShouldClose()) {
        float dt = GetFrameTime();
        draw(client, game);
        handle_controls(client, game);
        if (game->frame % 60 == 0) {
            for (int i = 1; i < game->n_entities; i++) {
                Vector2 target = (Vector2) {rand() / (float)RAND_MAX - 0.5,
                                            rand() / (float)RAND_MAX - 0.5};
                target = Vector2Normalize(target);
                game->entities[i].target_velocity = target;
            }
        }
        process_collisions(game, dt);
        if (game->elapsed_time - game->last_regenned > REGEN_EVERY_S) {
            regen_board(game);
            game->last_regenned = game->elapsed_time;
        }
        for (int i = game->n_entities - 1; i >= 0; i--) {
            // game->entities[i].size -= SIZE_DECAY_PER_S * dt;
            if (game->entities[i].size >= DEATH_SIZE) {
                if (!game->entities[i].attracted) {
                    game->entities[i].size = exp2(log2(game->entities[i].size) -
                                                  LOG_SIZE_DECAY_PER_S * dt);
                } else {
                    game->entities[i].size =
                        exp2(log2(game->entities[i].size) -
                             LOG_ATTRACTION_DECAY_PER_S * dt);
                }
            }
            if (game->entities[i].size < DEATH_SIZE) {
                destroy_entity(game, &game->entities[i]);
            }
        }
        generate_wall(client, game);
        game->frame++;
    }

    CloseWindow();
    free(client);
    return 0;
}

inline bool process_collisions(Game *game, float dt) {
    bool wall_changed = false;
    for (int i = 0; i < game->n_entities; i++) {
        Vector2 target = Vector2Scale(game->entities[i].target_velocity,
                                      MOVE_SPEED * (game->entities[i].attracted
                                          ? ATTRACTION_MULTIPLIER_VELOCITY
                                          : 1.0));
        game->entities[i].velocity = Vector2Add(
            game->entities[i].velocity,
            Vector2ClampValue(Vector2Subtract(target,
                                              game->entities[i].velocity),
                              0.0, ACCELERATION / deceleration(&game->entities[i]) * (game->entities[i].attracted
                                          ? ATTRACTION_MULTIPLIER_ACCELERATION
                                          : 1.0) * dt));
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
            if (!game->entities[i].attracted) {
                Vector2 epicenter = game->entities[i].position;
                int destroyed = 0;
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
                            destroyed++;
                        }
                    }
                }
                game->entities[i].size += destroyed * PER_DESTROYED_SIZE;
            } else {
                game->entities[i].size -= LOSE_BOUNCE_BONUS;
            }
        }
    }
    // for (int i = 0; i < 1; i++) {
    for (int i = 0; i < game->n_entities - 1; i++) {
        for (int j = i + 1; j < game->n_entities; j++) {
            if (Vector2Distance(game->entities[i].position, game->entities[j].position) >
                ent_radius(&game->entities[i]) + ent_radius(&game->entities[j])) {
                    continue;
            }
            float i_size = game->entities[i].size;
            float j_size = game->entities[j].size;
            if (i_size < DEATH_SIZE || j_size < DEATH_SIZE) {
                continue;
            }
            Vector2 direction = Vector2Normalize(Vector2Subtract(game->entities[i].position, game->entities[j].position));
            Vector2 i_velocity = game->entities[i].velocity;
            Vector2 j_velocity = game->entities[j].velocity;
            float i_dot = Vector2DotProduct(i_velocity, direction);
            float j_dot = Vector2DotProduct(j_velocity, direction);
            float to_redistribute = fmin(0.15, fmax(i_size, j_size));
            // split
            float i_kinetic_energy = i_size * i_dot * i_dot;
            float j_kinetic_energy = j_size * j_dot * j_dot;
            if (i_kinetic_energy + j_kinetic_energy == 0) {
                continue;
            }
            float i_fraction = 0.5 * (i_kinetic_energy - j_kinetic_energy) / fabsf(i_kinetic_energy - j_kinetic_energy);
            float j_fraction = 0.5 * (j_kinetic_energy - j_kinetic_energy) / fabsf(i_kinetic_energy - j_kinetic_energy);
            if (i_dot > 0 && j_dot > 0) {
                // both point in the j->i direction, therefore we give all the energy to j
                i_fraction = -0.5;
                j_fraction = 0.5;
            } else if (i_dot < 0 && j_dot < 0) {
                i_fraction = 0.5;
                j_fraction = -0.5;
            }
            // if (i == 0) {
            //     printf("to_redistribute: %f, i_fraction: %f, j_fraction: %f\n", to_redistribute, i_fraction, j_fraction);
            // }
            if (game->elapsed_time - game->entities[j].immune_since >= FIGHT_PERIOD
                && game->elapsed_time - game->entities[i].immune_since >= FIGHT_PERIOD) {
                game->entities[i].size += to_redistribute * i_fraction;
                game->entities[j].size += to_redistribute * j_fraction;
                game->entities[i].immune_since = game->elapsed_time;
                game->entities[j].immune_since = game->elapsed_time;
            }
            // add impulse
            i_velocity = Vector2Subtract(i_velocity, Vector2Scale(direction, i_dot));
            j_velocity = Vector2Subtract(j_velocity, Vector2Scale(direction, j_dot));
            i_velocity = Vector2Add(i_velocity, Vector2Scale(direction, sign(j_dot) * BALL_TO_BALL_IMPULSE));
            j_velocity = Vector2Add(j_velocity, Vector2Scale(direction, sign(i_dot) * BALL_TO_BALL_IMPULSE));
            game->entities[i].velocity = i_velocity;
            game->entities[j].velocity = j_velocity;
            game->entities[i].position = Vector2Add(game->entities[i].position, Vector2Scale(i_velocity, dt * 2));
            game->entities[j].position = Vector2Add(game->entities[j].position, Vector2Scale(j_velocity, dt * 2));
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
        new_camera.x = Clamp(new_camera.x, PI * 0.05, PI * 0.4);
        client->camera_angle = new_camera;
        client->camera_moved = true;
    }

    game->entities[0].target_velocity = (Vector2){0, 0};
    game->entities[0].attracted = false;

    Vector2 camera_forward = Vector2Normalize(Vector2Subtract(
        (Vector2){client->camera.target.x, client->camera.target.z},
        (Vector2){client->camera.position.x, client->camera.position.z}));
    // camera_forward = Vector2Scale(camera_forward, MOVE_SPEED);
    Vector2 camera_right = Vector2Rotate(camera_forward, PI / 2);
    Vector2 direction = Vector2Zero();
    if (IsKeyDown(KEY_K)) {
        direction =
            Vector2Add(direction, camera_forward);
    }
    if (IsKeyDown(KEY_J)) {
        direction =
            Vector2Subtract(direction, camera_forward);
    }
    if (IsKeyDown(KEY_L)) {
        direction =
            Vector2Add(direction, camera_right);
    }
    if (IsKeyDown(KEY_H)) {
        direction =
            Vector2Subtract(direction, camera_right);
    }
    if (Vector2Length(direction) > 0.0) {
        game->entities[0].target_velocity = Vector2Normalize(direction);
    } else {
        game->entities[0].target_velocity = (Vector2){0, 0};
    }
    if (IsKeyDown(KEY_SLASH)) {
        attract_to_nearest(&game->entities[0], game);
    }
}

void draw(Client *client, Game *game) {
    BeginDrawing();

    ClearBackground(BLACK);

    BeginMode3D(client->camera);

    BeginShaderMode(client->shader);

    Material sky_mat = LoadMaterialDefault();
    sky_mat.shader = client->shader;
    sky_mat.maps[MATERIAL_MAP_DIFFUSE].color = (Color){0, 0, 255, 255};
    DrawMesh(client->sky_wall, sky_mat, MatrixMultiply(
        MatrixRotateX(PI / 2), MatrixTranslate(0, SKY_SIZE / 2, -SKY_SIZE / 2)));
    DrawMesh(client->sky_wall, sky_mat, MatrixMultiply(
        MatrixRotateX(-PI / 2), MatrixTranslate(0, SKY_SIZE / 2, SKY_SIZE / 2)));
    DrawMesh(client->sky_wall, sky_mat,
             MatrixMultiply(
                 MatrixMultiply(MatrixRotateX(PI / 2), MatrixRotateY(-PI / 2)),
                 MatrixTranslate(SKY_SIZE / 2, SKY_SIZE / 2, 0)));
    DrawMesh(client->sky_wall, sky_mat,
             MatrixMultiply(
                 MatrixMultiply(MatrixRotateX(PI / 2), MatrixRotateY(PI / 2)),
                 MatrixTranslate(-SKY_SIZE / 2, SKY_SIZE / 2, 0)));
    DrawMesh(client->sky_wall, sky_mat,
            MatrixMultiply(MatrixRotateX(PI), MatrixTranslate(0, SKY_SIZE / 2, 0)));

    DrawModel(client->outer_walls, (Vector3){0.0f, 0.0f, 0.0f}, 1.0f,
              (Color){255, 255, 0, 255});

    DrawPlane((Vector3){0.0f, 0.0f, 0.0f}, (Vector2){GROUND_EXTENTS, GROUND_EXTENTS}, (Color){0, 255, 0, 255});

    Material mat = LoadMaterialDefault();
    mat.maps[MATERIAL_MAP_DIFFUSE].color = (Color){0, 255, 255, 0};
    for (int y = 0; y < BOARD_HEIGHT; y += RENDER_CHUNK) {
        for (int x = 0; x < BOARD_WIDTH; x += RENDER_CHUNK) {
            if (client->walls[y / RENDER_CHUNK][x / RENDER_CHUNK].meshCount == 0) {
                continue;
            }
            DrawModel(client->walls[y / RENDER_CHUNK][x / RENDER_CHUNK],
                      (Vector3){x - BOARD_WIDTH / 2, 0.0f, y - BOARD_HEIGHT / 2},
                      1.0f, WHITE);
        }
    }
    // DrawModel(client->walls, (Vector3){-BOARD_WIDTH / 2 + 1, 0.0f, -BOARD_HEIGHT / 2 + 1}, 1.0f, WHITE);

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
    // snprintf(textbuf, 256,
    //          "x: %.1f, y: %.1f, vx: %.1f, vy: %.1f, tvx: %.1f, tvy: %.1f",
    //          game->entities[0].position.x, game->entities[0].position.y,
    //          game->entities[0].velocity.x, game->entities[0].velocity.y,
    //          game->entities[0].target_velocity.x,
    //          game->entities[0].target_velocity.y);
    snprintf(textbuf, 256,
        "size: %.2f", game->entities[0].size);
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
        float radius = ent_radius(&game->entities[i]);
        int radius_ceil = ceil(radius * 1.2);
        for (int oy = -radius_ceil; oy <= radius_ceil; oy++) {
            for (int ox = -radius_ceil; ox <= radius_ceil; ox++) {
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
    unsigned char wall_data[RENDER_CHUNK][RENDER_CHUNK][3] = {0};
    Image walls = (Image){.width = RENDER_CHUNK,
                          .height = RENDER_CHUNK,
                          .format = PIXELFORMAT_UNCOMPRESSED_R8G8B8,
                          .data = (void *)wall_data};
    // int have_changed = 0;
    for (int y = 0; y < BOARD_HEIGHT; y += RENDER_CHUNK) {
        for (int x = 0; x < BOARD_WIDTH; x += RENDER_CHUNK) {
            bool has_changed = false;
            for (int oy = 0; oy < RENDER_CHUNK; oy++) {
                for (int ox = 0; ox < RENDER_CHUNK; ox++) {
                    if (game->board[y + oy][x + ox] != client->past_board[y + oy][x + ox]) {
                        has_changed = true;
                        break;
                    }
                }
                if (has_changed) {
                    break;
                }
            }
            if (!has_changed) {
                continue;
            }
            if (client->walls[y / RENDER_CHUNK][x / RENDER_CHUNK].meshCount > 0) {
                UnloadModel(client->walls[y / RENDER_CHUNK][x / RENDER_CHUNK]);
            }
            // have_changed++;
            for (int oy = 0; oy < RENDER_CHUNK; oy++) {
                for (int ox = 0; ox < RENDER_CHUNK; ox++) {
                    unsigned char col = game->board[y + oy][x + ox] == WALL ? 255 : 0;
                    wall_data[oy][ox][0] = col;
                    wall_data[oy][ox][1] = col;
                    wall_data[oy][ox][2] = col;
                }
            }
            Mesh walls_mesh = GenMeshCubicmap(walls, (Vector3){1.0, WALL_HEIGHT, 1.0});
            client->walls[y / RENDER_CHUNK][x / RENDER_CHUNK] = LoadModelFromMesh(walls_mesh);
            client->walls[y / RENDER_CHUNK][x / RENDER_CHUNK].materials[0].shader = client->shader;
            client->walls[y / RENDER_CHUNK][x / RENDER_CHUNK].materials[0].maps[MATERIAL_MAP_DIFFUSE].color =
                (Color){0, 255, 255, 255};
        }
    }
    // float changed_frac = have_changed / (float)(BOARD_HEIGHT / RENDER_CHUNK * BOARD_WIDTH / RENDER_CHUNK);
    // printf("changed_frac: %f\n", changed_frac);
    copy_past_board(client, game);
}

void copy_past_board(Client *client, Game *game) {
    for (int y = 0; y < BOARD_HEIGHT; y++) {
        for (int x = 0; x < BOARD_WIDTH; x++) {
            client->past_board[y][x] = game->board[y][x];
        }
    }
}

void destroy_entity(Game *game, Entity *e) {
    int offset = ((int)(e - game->entities));
    if (offset == 0) {
        printf("tried to destroy player\n");
        printf("size: %f\n", e->size);
        exit(1);
    }
    Entity tmp = game->entities[offset];
    game->entities[offset] = game->entities[game->n_entities - 1];
    game->entities[game->n_entities - 1] = tmp;
    game->n_entities--;
}

void attract_to_nearest(Entity *entity, Game *game) {
    if (entity->size < ATTRACT_SIZE) {
        return;
    }
    float nearest_distance = __FLT_MAX__;
    int nearest_index = -1;
    for (int i = 0; i < game->n_entities; i++) {
        if (entity == &game->entities[i]) {
            continue;
        }
        if (game->entities[i].size < DEATH_SIZE) {
            continue;
        }
        if (Vector2Distance(entity->position, game->entities[i].position) < ent_radius(entity) + ent_radius(&game->entities[i])) {
            continue;
        }
        float distance = Vector2Distance(entity->position, game->entities[i].position);
        if (distance < nearest_distance) {
            nearest_distance = distance;
            nearest_index = i;
        }
    }
    Vector2 direction = Vector2Normalize(Vector2Subtract(game->entities[nearest_index].position, entity->position));
    entity->attracted = true;
    entity->target_velocity = direction;
}
