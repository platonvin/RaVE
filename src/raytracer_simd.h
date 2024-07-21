#pragma once
#ifndef SIMD_RAVE_H
#define SIMD_RAVE_H

#include <assert.h>
// #include "../dependencies/sl_vec.h"
// #include <math.h>

// float _modf( float arg ) {
//     float integral_part;
//     return modff(arg, &integral_part);
// }
#include <stdbool.h>
#define println printf("%s:%d %s\n",__FILE__, __LINE__, __FUNCTION__);

#include <immintrin.h>
__m256 a;
// #define fract(x) _modf(x);
#define PI 3.14159265358979323846
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct rave_vec2{float x,y    ;}rave_vec2;
typedef struct rave_vec3{float x,y,z  ;}rave_vec3;
typedef struct rave_vec4{float x,y,z,w;}rave_vec4;

typedef struct rave_ivec2{int x,y    ;}rave_ivec2;
typedef struct rave_ivec3{int x,y,z  ;}rave_ivec3;
typedef struct rave_ivec4{int x,y,z,w;}rave_ivec4;

typedef union rave_pack{
    __m256 body;
    struct {
        rave_vec4 v1, v2;
    } repr;
}rave_pack;

typedef struct simd_rave_ivec4{int x,y,z,w;}simd_rave_ivec4;

typedef struct Material {
    rave_vec3 color;
    float roughness;
    float emmitance;
} Material;

// typedef struct rave_color {
//     ravec3 rgb;
// }rave_color;
typedef uint8_t voxel; //TODO

typedef struct rave_ctx {
    voxel  (*rave_get_voxel)(int x, int y, int z);
    rave_vec3 (*rave_get_ray_pos)(int local_x, int local_y, int local_z);
    rave_vec3 (*rave_get_ray_dir)(int local_x, int local_y, int local_z);
    void   (*rave_store_light)(int local_x, int local_y, int local_z, rave_vec3 light);
    Material (*rave_get_material)(voxel voxel);
    uint32_t max_steps;
    uint32_t max_reflections;
    uint32_t bounds;
}rave_ctx; 

float dot(rave_vec3 v1, rave_vec3 v2);

rave_vec3 reflect(rave_vec3 I, rave_vec3 N);
rave_vec3 normalize(rave_vec3 v);
rave_vec3 mix(rave_vec3 v1, rave_vec3 v2, float t);
float distance(rave_vec3 v1, rave_vec3 v2);
rave_vec3 add(rave_vec3 v1, rave_vec3 v2);
rave_vec3 sub(rave_vec3 v1, rave_vec3 v2);
rave_vec3 mul(rave_vec3 v, float scalar);
rave_vec3 clamp(rave_vec3 v, float minVal, float maxVal);
rave_vec3 cross(rave_vec3 v1, rave_vec3 v2);
void rave_dispatch(rave_ctx* ctx, int count_x, int count_y, int count_z);
void rave_sync();

void pf(rave_vec3 v);
void pi(rave_ivec3 v);

#ifdef __cplusplus
}
#endif

#endif // SIMD_RAVE_H