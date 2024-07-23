#pragma once
#ifndef SIMD_RAVE_H
#define SIMD_RAVE_H

#include <assert.h>
#include <stdint.h>
#include <stdbool.h>

#define println printf("%s:%d %s\n",__FILE__, __LINE__, __FUNCTION__);
#define PI 3.14159265358979323846

#ifdef __cplusplus
extern "C" {
#endif

// typedef struct rave_vec2{float x,y    ;}rave_vec2;
typedef struct rave_vec3{float x,y,z  ;}rave_vec3;
typedef struct rave_vec4{
    float x,y,z,w;
}rave_vec4;

// typedef struct rave_ivec2{int x,y    ;}rave_ivec2;
typedef struct rave_ivec3{int x,y,z  ;}rave_ivec3;
// typedef struct rave_ivec4{
//     // struct {int vec[4];};
//     struct {int x,y,z,w;};
// }rave_ivec4;

typedef struct Material {
    rave_vec3 color;
    float roughness;
    float emmitance;
} Material;

typedef uint8_t voxel; //TODO

// #ifndef

voxel  rave_get_voxel(int x, int y, int z);
rave_vec3 rave_get_ray_pos(int local_x, int local_y, int local_z);
rave_vec3 rave_get_ray_dir(int local_x, int local_y, int local_z);
void   rave_store_light(int local_x, int local_y, int local_z, rave_vec3 light);
Material rave_get_material(voxel voxel);
extern uint32_t max_steps;
extern uint32_t max_reflections;
// extern uint32_t bounds;

typedef struct rave_ctx {
    // voxel  (*rave_get_voxel)(int x, int y, int z);
    // rave_vec3 (*rave_get_ray_pos)(int local_x, int local_y, int local_z);
    // rave_vec3 (*rave_get_ray_dir)(int local_x, int local_y, int local_z);
    // void   (*rave_store_light)(int local_x, int local_y, int local_z, rave_vec3 light);
    // Material (*rave_get_material)(voxel voxel);
    // uint32_t max_steps;
    // uint32_t max_reflections;
    // uint32_t bounds;
}rave_ctx; 

float dot3(rave_vec3 v1, rave_vec3 v2);
rave_vec3 reflect3(rave_vec3 I, rave_vec3 N);
rave_vec3 normalize3(rave_vec3 v);
rave_vec3 mix3(rave_vec3 v1, rave_vec3 v2, float t);
float distance3(rave_vec3 v1, rave_vec3 v2);
rave_vec3 add3(rave_vec3 v1, rave_vec3 v2);
rave_vec3 sub3(rave_vec3 v1, rave_vec3 v2);
rave_vec3 mul3(rave_vec3 v, float scalar);
rave_vec3 clamp3(rave_vec3 v, float minVal, float maxVal);
rave_vec3 cross3(rave_vec3 v1, rave_vec3 v2);


void rave_dispatch_simd(rave_ctx* ctx, int count_x, int count_y, int count_z);
void rave_sync_simd();
void rave_init();

#define pf(v) {\
float* ptr = &v;\
printf(#v "= %f %f %f %f\n", ptr[0], ptr[1], ptr[2], ptr[3]);\
}
#define pi(v) {\
int* ptr = &v;\
printf(#v "= %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3]);\
}

#ifdef __cplusplus
}
#endif

#endif // SIMD_RAVE_H