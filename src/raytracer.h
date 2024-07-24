#pragma once
#ifndef RAVE_H
#define RAVE_H

/*
to use __custom__type__ instead of uint8_t 
#define RAVE_CUSTOM_VOXEL_TYPE __custom__type__ 
#define RAVE_EMPTY_VOXEL __custom__type__empty__voxel__value__

to disable SIMD 
#define RAVE_NO_SIMD
*/
#include <assert.h>
#include <stdint.h>
#include <stdbool.h>

#define println printf("%s:%d %s\n",__FILE__, __LINE__, __FUNCTION__);
#define RAVE_PI 3.14159265358979323846

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

#ifdef RAVE_CUSTOM_VOXEL_TYPE
typedef RAVE_CUSTOM_VOXEL_TYPE rave_voxel;
#else
typedef uint8_t rave_voxel;
#endif

#ifdef RAVE_CUSTOM_VOXEL_TYPE
// rave_empty_voxel = RAVE_EMPTY_VOXEL;
#else
// rave_empty_voxel = 0;
#endif

rave_voxel rave_get_voxel(int x, int y, int z);
rave_vec3  rave_get_ray_pos(int local_x, int local_y, int local_z);
rave_vec3  rave_get_ray_dir(int local_x, int local_y, int local_z);
void       rave_store_light(int local_x, int local_y, int local_z, rave_vec3 light);
Material   rave_get_material(rave_voxel voxel);
extern uint32_t max_steps;
extern uint32_t max_reflections;
extern rave_voxel rave_empty_voxel;

float     rave_dot3(rave_vec3 v1, rave_vec3 v2);
rave_vec3 rave_reflect3(rave_vec3 I, rave_vec3 N);
rave_vec3 rave_normalize3(rave_vec3 v);
rave_vec3 rave_mix3(rave_vec3 v1, rave_vec3 v2, float t);
float     rave_distance3(rave_vec3 v1, rave_vec3 v2);
rave_vec3 rave_add3(rave_vec3 v1, rave_vec3 v2);
rave_vec3 rave_sub3(rave_vec3 v1, rave_vec3 v2);
rave_vec3 rave_mul3(rave_vec3 v, float scalar);
rave_vec3 rave_clamp3(rave_vec3 v, float minVal, float maxVal);
rave_vec3 rave_cross3(rave_vec3 v1, rave_vec3 v2);

void rave_dispatch(int count_x, int count_y, int count_z);
void rave_init();

//TODO: reuse pthreads, explicit sync
// void rave_sync();

//printf __m128
#define __pf(v) {\
float* ptr = &v;\
printf(#v "= %f %f %f %f\n", ptr[0], ptr[1], ptr[2], ptr[3]);\
}
#define __ti(v) {\
int* ptr = &v;\
printf(#v "= %d %d %d %d\n", ptr[0], ptr[1], ptr[2], ptr[3]);\
}

#ifdef __cplusplus
}
#endif

#endif // RAVE_H