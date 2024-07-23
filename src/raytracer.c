#include "raytracer.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <immintrin.h>

// static unsigned long int next = 1;


float Random1D(){
    float res = ((float)rand()) / ((float)RAND_MAX);
        //   res = 2.0*sinf(res*523145.1245) - 1.0;
          res = 2.0*res - 1.0;
    return (res);
}

rave_vec4 Random4D() {
    float a,b,c;
    a = Random1D();
    b = Random1D();
    c = Random1D();
    return (rave_vec4){a,b,c, 0};
}

rave_vec4 random_sphere_point_simd(rave_vec4 rand) {
    float ang1 = (rand.x + 1.0) * PI; // [-1..1) -> [0..2*PI)
    float u = rand.y; // [-1..1), cos and acos(2v-1) cancel each other out, so we arrive at [-1..1)
    float u2 = u * u;
    float sqrt1MinusU2 = sqrtf(1.0 - u2);
    float x = sqrt1MinusU2 * cosf(ang1);
    float y = sqrt1MinusU2 * sinf(ang1);
    float z = u;
    return (rave_vec4){x, y, z, 0};
}

rave_vec4 normal_oriented_hemisphere_point_simd(rave_vec4 rand, rave_vec4 n){
    rave_vec4 v = random_sphere_point_simd(rand);
    float dot = 
          v.x * n.x +
          v.y * n.y +
          v.z * n.z;
    float sign = (dot>0) ? 1.0 : -1.0;
    v.x *= sign;
    v.y *= sign;
    v.z *= sign;
    
    //for simd
    // v.w *= sign;
    // v.vec[3] *= sign;
    return v;
}
inline __m128 _mm_abs_ps(__m128 m) {
    return _mm_andnot_ps(_mm_set1_ps(-0.0f), m);
}
inline bool initTvals_simd(__m128* restrict tMax, __m128* restrict tDelta, __m128i* restrict blockPos, rave_vec4 rayOrigin, rave_vec4 rayDirection){


    // rave_vec4 effective_origin = rayOrigin;
    __m128 _eo = _mm_set_ps(
            rayOrigin.w,
            rayOrigin.z,
            rayOrigin.y,
            rayOrigin.x
        );
    // pf(_eo);
        // block_corner1.x = (floorf(effective_origin.x) - effective_origin.x)/rayDirection.x; //now corners are relative vectors
        // block_corner1.y = (floorf(effective_origin.y) - effective_origin.y)/rayDirection.y; //now corners are relative vectors
        // block_corner1.z = (floorf(effective_origin.z) - effective_origin.z)/rayDirection.z; //now corners are relative vectors
        // block_corner1.w = (floorf(effective_origin.w) - effective_origin.w)/rayDirection.w; //now corners are relative vectors
    __m128 _rd = _mm_set_ps(
            rayDirection.w,
            rayDirection.z,
            rayDirection.y,
            rayDirection.x
        );
    __m128 block_corner1 = 
        _mm_div_ps(
            _mm_sub_ps(
                _mm_floor_ps(_eo),
                _eo
            ),
            _rd
        );
    // pf(block_corner1);

    //   block_corner2.x = (floor(effective_origin.x) - effective_origin.x)/rayDirection.x  + 1.0/rayDirection.x;
    //   block_corner2.y = (floor(effective_origin.y) - effective_origin.y)/rayDirection.y  + 1.0/rayDirection.y;
    //   block_corner2.z = (floor(effective_origin.z) - effective_origin.z)/rayDirection.z  + 1.0/rayDirection.z;
    //   block_corner2.w = (floor(effective_origin.w) - effective_origin.w)/rayDirection.w  + 1.0/rayDirection.w;
    __m128 block_corner2 = 
        _mm_add_ps(
            block_corner1,
            _mm_div_ps(
                _mm_set1_ps(1.0),
                _rd
            )
        );

    // (*tMax).x = fmaxf(block_corner1.x, block_corner2.x); //1 of theese will be negative so max is just to get positive
    // (*tMax).y = fmaxf(block_corner1.y, block_corner2.y);
    // (*tMax).z = fmaxf(block_corner1.z, block_corner2.z);
    // (*tMax).w = fmaxf(block_corner1.w, block_corner2.w);
    (*tMax) = _mm_max_ps(block_corner1, block_corner2);

    // (*tDelta).x = 1.0 / fabsf(rayDirection.x); //how many dir vectors needeed to move 1.0 across each axys
    // (*tDelta).y = 1.0 / fabsf(rayDirection.y); //how many dir vectors needeed to move 1.0 across each axys
    // (*tDelta).z = 1.0 / fabsf(rayDirection.z); //how many dir vectors needeed to move 1.0 across each axys
    // (*tDelta).w = 1.0 / fabsf(rayDirection.w); //how many dir vectors needeed to move 1.0 across each axys
    (*tDelta) = _mm_div_ps(
        _mm_set1_ps(1.0),
        _mm_abs_ps(_rd));

    // _MM_FROUND_TRUNC should be used instead of _MM_FROUND_FLOOR
    (*blockPos) = _mm_cvtps_epi32(_mm_round_ps((_eo), (_MM_FROUND_RAISE_EXC | _MM_FROUND_TO_POS_INF)));
    // (*blockPos) = (rave_vec4){effective_origin.x, effective_origin.y, effective_origin.z, effective_origin.w}; //round origin to block pos

    return true;
}
bool initTvals(rave_vec3* tMax, rave_vec3* tDelta, rave_ivec3* blockPos, rave_vec3 rayOrigin, rave_vec3 rayDirection){
    rave_vec3 effective_origin = rayOrigin;
    // printf("effective_origin= %f %f %f\n", effective_origin.x, effective_origin.y, effective_origin.z);

    rave_vec3 block_corner1;
              block_corner1.x = (floorf(effective_origin.x) - effective_origin.x)/rayDirection.x; //now corners are relative vectors
              block_corner1.y = (floorf(effective_origin.y) - effective_origin.y)/rayDirection.y; //now corners are relative vectors
              block_corner1.z = (floorf(effective_origin.z) - effective_origin.z)/rayDirection.z; //now corners are relative vectors
    // printf("block_corner1/= %f %f %f\n", block_corner1.x, block_corner1.y, block_corner1.z);

    rave_vec3 block_corner2;
              block_corner2.x = (floorf(effective_origin.x) - effective_origin.x)/rayDirection.x  + 1.0/rayDirection.x;
              block_corner2.y = (floorf(effective_origin.y) - effective_origin.y)/rayDirection.y  + 1.0/rayDirection.y;
              block_corner2.z = (floorf(effective_origin.z) - effective_origin.z)/rayDirection.z  + 1.0/rayDirection.z;

    (*tMax).x = fmaxf(block_corner1.x, block_corner2.x); //1 of theese will be negative so max is just to get positive
    (*tMax).y = fmaxf(block_corner1.y, block_corner2.y);
    (*tMax).z = fmaxf(block_corner1.z, block_corner2.z);

    (*tDelta).x = 1.0 / fabsf(rayDirection.x); //how many dir vectors needeed to move 1.0 across each axys
    (*tDelta).y = 1.0 / fabsf(rayDirection.y); //how many dir vectors needeed to move 1.0 across each axys
    (*tDelta).z = 1.0 / fabsf(rayDirection.z); //how many dir vectors needeed to move 1.0 across each axys

    (*blockPos) = (rave_ivec3){ceilf(effective_origin.x), ceilf(effective_origin.y), ceilf(effective_origin.z)}; //round origin to block pos

    return true;
}
static __m128i _my_and_epi32(__m128i __a, __m128i __b)
{
  return (__m128i)((__v4su)__a & (__v4su)__b);
}

int CastRay_precise(rave_ctx* ctx, rave_vec3 rayOrigin, rave_vec3 rayDirection, 
        float* fraction, rave_vec3* normal, Material* material){
    bool block_hit = false;

    rave_ivec3 steps;
    
    steps = (rave_ivec3){(rayDirection.x > 0), (rayDirection.y > 0), (rayDirection.z > 0)};
    steps.x = 2 * steps.x - 1;
    steps.y = 2 * steps.y - 1;
    steps.z = 2 * steps.z - 1;

    rave_vec3 fsteps = {steps.x, steps.y, steps.z};

    rave_vec3 tMax = {0};
    rave_vec3 tDelta = {0};
    rave_ivec3 voxel_pos = {0};
    float block_fraction = 0.0;

    rave_ivec3 fcurrentStepDiretion = {0};

    initTvals(&tMax, &tDelta, &voxel_pos, rayOrigin, rayDirection); //does not intersect with scene
    // printf("tMax= %f %f %f\n", tMax.x, tMax.y, tMax.z);
    // printf("tDelta= %f %f %f\n", tDelta.x, tDelta.y, tDelta.z);
    // printf("voxel_pos= %d %d %d\n", voxel_pos.x, voxel_pos.y, voxel_pos.z);
    
    int current_voxel = rave_get_voxel(voxel_pos.x, voxel_pos.y, voxel_pos.z);

    assert(max_steps != 0);

    int iterations = 0;

    while (true) {
        fcurrentStepDiretion = (rave_ivec3){0};
        if(tMax.x <= tMax.y){
            if(tMax.z <= tMax.x){
                fcurrentStepDiretion.z = 1;
                voxel_pos.z += steps.z;
                tMax.z += tDelta.z;
            } else {
                fcurrentStepDiretion.x = 1;        
                voxel_pos.x += steps.x;
                tMax.x += tDelta.x;
            }
        } else{
            if (tMax.y <= tMax.z) {
                fcurrentStepDiretion.y = 1;        
                voxel_pos.y += steps.y;
                tMax.y += tDelta.y;
            } else {
                fcurrentStepDiretion.z = 1;        
                voxel_pos.z += steps.z;
                tMax.z += tDelta.z;
            }
        }
        // bool xLy = tMax.x <= tMax.y;
        // bool xLz = tMax.x <= tMax.z;
        // bool yLz = tMax.y <= tMax.z;

        // fcurrentStepDiretion.x = (int)((int)(( xLy) && ( xLz)));
        // fcurrentStepDiretion.y = (int)((int)((!xLy) && ( yLz)));
        // fcurrentStepDiretion.z = (int)((int)((!xLz) && (!yLz)));
        
        // voxel_pos.x += steps.x * fcurrentStepDiretion.x;
        // voxel_pos.y += steps.y * fcurrentStepDiretion.y;
        // voxel_pos.z += steps.z * fcurrentStepDiretion.z;
        
        // tMax.x += tDelta.x * fcurrentStepDiretion.x;
        // tMax.y += tDelta.y * fcurrentStepDiretion.y;
        // tMax.z += tDelta.z * fcurrentStepDiretion.z;
        
        current_voxel = rave_get_voxel(voxel_pos.x, voxel_pos.y, voxel_pos.z);
        
        if (current_voxel != 0){
            block_hit = true;
            // println
            break;
        }
        if ((iterations++ >= max_steps)) {
            block_hit = false;
            break;
        }
    }


    (*normal).x = -(fsteps.x * fcurrentStepDiretion.x);
    (*normal).y = -(fsteps.y * fcurrentStepDiretion.y);
    (*normal).z = -(fsteps.z * fcurrentStepDiretion.z);
    // printf("(*normal)= %f %f %f\n", (*normal).x, (*normal).y, (*normal).z);

    rave_vec3 tFinal;
              tFinal.x = tMax.x - tDelta.x;
              tFinal.y = tMax.y - tDelta.y;
              tFinal.z = tMax.z - tDelta.z;
    // block_fraction += dot(tFinal, fcurrentStepDiretion);
    // printf("tMax= %f %f %f\n", tMax.x, tMax.y, tMax.z);
    // printf("tDelta= %f %f %f\n", tDelta.x, tDelta.y, tDelta.z);
    // printf("tFinal= %f %f %f\n", tFinal.x, tFinal.y, tFinal.z);
    
    block_fraction = 
        tFinal.x * (float)fcurrentStepDiretion.x+
        tFinal.y * (float)fcurrentStepDiretion.y+
        tFinal.z * (float)fcurrentStepDiretion.z;

    // printf("voxel %d\n", current_voxel);
    (*material) = rave_get_material(current_voxel);
    (*fraction) = block_fraction;

    return (iterations);
}

float dot3(rave_vec3 v1, rave_vec3 v2) {
    return v1.x*v2.x + 
           v1.y*v2.y + 
           v1.z*v2.z;
}
rave_vec3 reflect3(rave_vec3 I, rave_vec3 N) {
    float dotProduct = dot3(I, N);
    rave_vec3 result;
    result.x = I.x - 2.0f * dotProduct * N.x;
    result.y = I.y - 2.0f * dotProduct * N.y;
    result.z = I.z - 2.0f * dotProduct * N.z;
    return result;
}
rave_vec3 normalize3(rave_vec3 v) {
    float length = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
    rave_vec3 result;
    if (length > 0.0f) {
        result.x = v.x / length;
        result.y = v.y / length;
        result.z = v.z / length;
    } else {
        result.x = result.y = result.z = 0.0f;
    }
    return result;
}
rave_vec3 cross3(rave_vec3 v1, rave_vec3 v2) {
    rave_vec3 result;
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}
rave_vec3 mix3(rave_vec3 v1, rave_vec3 v2, float t) {
    rave_vec3 result;
    result.x = v1.x * (1.0f - t) + v2.x * t;
    result.y = v1.y * (1.0f - t) + v2.y * t;
    result.z = v1.z * (1.0f - t) + v2.z * t;
    return result;
}
float distance3(rave_vec3 v1, rave_vec3 v2) {
    float dx = v1.x - v2.x;
    float dy = v1.y - v2.y;
    float dz = v1.z - v2.z;
    return sqrtf(dx * dx + dy * dy + dz * dz);
}
rave_vec3 add3(rave_vec3 v1, rave_vec3 v2) {
    rave_vec3 result;
    result.x = v1.x + v2.x;
    result.y = v1.y + v2.y;
    result.z = v1.z + v2.z;
    return result;
}
rave_vec3 sub3(rave_vec3 v1, rave_vec3 v2) {
    rave_vec3 result;
    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    result.z = v1.z - v2.z;
    return result;
}
rave_vec3 mul3(rave_vec3 v, float scalar) {
    rave_vec3 result;
    result.x = v.x * scalar;
    result.y = v.y * scalar;
    result.z = v.z * scalar;
    return result;
}
rave_vec3 clamp3(rave_vec3 v, float minVal, float maxVal) {
    rave_vec3 result;
    result.x = fminf(fmaxf(v.x, minVal), maxVal);
    result.y = fminf(fmaxf(v.y, minVal), maxVal);
    result.z = fminf(fmaxf(v.z, minVal), maxVal);
    return result;
}

float dot(rave_vec4 v1, rave_vec4 v2) {
    // return v1.x*v2.x + 
    //        v1.y*v2.y + 
    //        v1.z*v2.z;
}
// rave_vec4 reflect_simd(rave_vec4 I, rave_vec4 N) {
//     float dotProduct = dot(I, N);
//     rave_vec4 result = {};
//     result.x = I.x - 2.0f * dotProduct * N.x;
//     result.y = I.y - 2.0f * dotProduct * N.y;
//     result.z = I.z - 2.0f * dotProduct * N.z;
//     return result;
// }
// rave_vec4 normalize(rave_vec4 v) {
//     float length = sqrtf(v.x*v.x + 
//                          v.y*v.y + 
//                          v.z*v.z);
                         
//     rave_vec4 result = {};

//     assert(length != 0.0f);
//     assert(length > 0.0f);
//     if (length > 0.0f) {
//         result.x = v.x / length;
//         result.y = v.y / length;
//         result.z = v.z / length;
//     } else {
//         // abort();
//         // result = (rave_vec4){};
//     }
//     return result;
// }
// rave_vec4 cross(rave_vec4 v1, rave_vec4 v2) {
//     rave_vec4 result;
//     result.x = v1.y * v2.z - v1.z * v2.y;
//     result.y = v1.z * v2.x - v1.x * v2.z;
//     result.z = v1.x * v2.y - v1.y * v2.x;
//     return result;
// }
// rave_vec4 mix(rave_vec4 v1, rave_vec4 v2, float t) {
//     rave_vec4 result;
//     assert(t >= 0.0);
//     assert(t <= 1.0);
//     result.x = v1.x * (1.0f - t) + v2.x * t;
//     result.y = v1.y * (1.0f - t) + v2.y * t;
//     result.z = v1.z * (1.0f - t) + v2.z * t;
//     // result.w = v1.w * (1.0f - t) + v2.w * t;
//     return result;
// }
// float distance(rave_vec4 v1, rave_vec4 v2) {
//     float dx = v1.x - v2.x;
//     float dy = v1.y - v2.y;
//     float dz = v1.z - v2.z;
//     // float dw = v1.w - v2.w;
//     return sqrtf(dx * dx + dy * dy + dz * dz);
// }
// rave_vec4 add(rave_vec4 v1, rave_vec4 v2) {
//     rave_vec4 result;
//     result.x = v1.x + v2.x;
//     result.y = v1.y + v2.y;
//     result.z = v1.z + v2.z;
//     result.w = v1.w + v2.w;
//     return result;
// }
// rave_vec4 sub(rave_vec4 v1, rave_vec4 v2) {
//     rave_vec4 result;
//     result.x = v1.x - v2.x;
//     result.y = v1.y - v2.y;
//     result.z = v1.z - v2.z;
//     result.w = v1.w - v2.w;
//     return result;
// }
// rave_vec4 mul(rave_vec4 v, float scalar) {
//     rave_vec4 result;
//     result.x = v.x * scalar;
//     result.y = v.y * scalar;
//     result.z = v.z * scalar;
//     result.w = v.w * scalar;
//     return result;
// }
// rave_vec4 clamp(rave_vec4 v, float minVal, float maxVal) {
//     rave_vec4 result;
//     result.x = fminf(fmaxf(v.x, minVal), maxVal);
//     result.y = fminf(fmaxf(v.y, minVal), maxVal);
//     result.z = fminf(fmaxf(v.z, minVal), maxVal);
//     result.w = fminf(fmaxf(v.w, minVal), maxVal);
//     return result;
// }
// void pf(rave_vec4 v) {
//     printf("rave_vec4: (%f, %f, %f)\n", v.x, v.y, v.z);
// }
// void pi(rave_ivec3 v) {
//     printf("rave_vec4: (%d, %d, %d)\n", v.x, v.y, v.z);
// }
void ProcessHit_simd(rave_vec4* origin, rave_vec4* direction, 
                float fraction, rave_vec4 normal, Material material, 
                rave_vec4* accumulated_light, rave_vec4* accumulated_reflection){

            rave_vec4 new_origin;
                      new_origin.x = (*origin).x + (fraction * (*direction).x);
                      new_origin.y = (*origin).y + (fraction * (*direction).y);
                      new_origin.z = (*origin).z + (fraction * (*direction).z);
// println
            normal = normalize(normal);

            rave_vec4 rnd = Random4D();
            
            rave_vec4 hemisphereDistributedDirection = normal_oriented_hemisphere_point_simd(rnd, normal);
// println 

            rave_vec4 new_direction = normalize(hemisphereDistributedDirection);
            // rave_vec4 new_direction = (hemisphereDistributedDirection);

            bool ass = fabsf(normal.x + normal.y + normal.z) == 1.0;
            if(!ass){
                printf("normal->x %f\n", normal.x);
                printf("normal->y %f\n", normal.y);
                printf("normal->z %f\n", normal.z);
                printf("\n");
                abort();
            }
            // assert();

            // rave_vec4 idealReflection = (reflect((*direction), normal));
// println
            rave_vec4 idealReflection = normalize(reflect_simd((*direction), normal));

// println
            new_direction = normalize(mix(idealReflection, new_direction, material.roughness));
            new_origin.x += normal.x * 0.01;
            new_origin.y += normal.y * 0.01;
            new_origin.z += normal.z * 0.01;

            (*accumulated_reflection).x *= material.color.x;
            (*accumulated_reflection).y *= material.color.y;
            (*accumulated_reflection).z *= material.color.z;
            
            (*accumulated_light).x += (*accumulated_reflection).x * material.emmitance;
            (*accumulated_light).y += (*accumulated_reflection).y * material.emmitance;
            (*accumulated_light).z += (*accumulated_reflection).z * material.emmitance;
            // accumulated_light += vec3(0.8) * 0.1 * accumulated_reflection;

            // direction = reflect(direction,normal);
            (*direction) = new_direction;
            (*origin) = new_origin;

}

static rave_vec4 trace_ray(rave_ctx* ctx, rave_vec4 origin, rave_vec4 direction){
    //finds intersections and processes hits
    float fraction = 0.0;
    rave_vec4 normal = {0};
    Material material = {0};

    rave_vec4 light = {0};
    rave_vec4 reflection = {1,1,1,0};
// println
    for (int i = 0; (i < max_reflections); i++){
        
        int hit=1;
        hit   = CastRay_precise(ctx, (rave_vec3){origin.x, origin.y, origin.z}, (rave_vec3){direction.x, direction.y, direction.z}, &fraction, &normal, &material);
        // printf("fraction= %f\n", fraction);
        // printf("material= %f %f %f _ %f %f\n", material.color.x, material.color.y, material.color.z, material.emmitance, material.roughness);
        // printf("normal= %f %f %f %f\n", normal.x, normal.y, normal.z, normal.w);
        // printf("hit= %d\n", hit);
        int hit_s=1;
        // hit_s = CastRay_precise_simd(ctx, origin, direction, &fraction, &normal, &material);
        // printf("SIMD fraction= %f\n", fraction);
        // printf("SIMD material= %f %f %f _ %f %f\n", material.color.x, material.color.y, material.color.z, material.emmitance, material.roughness);
        // printf("SIMD normal= %f %f %f %f\n", normal.x, normal.y, normal.z, normal.w);
        // printf("SIMD hit= %d\n", hit_s);
        // printf("\n");
        if(!hit) break;

        // if(!hit_s) break;

        ProcessHit_simd(&origin, &direction, 
            fraction, normal, material, 
            &light, &reflection);
        
        
        // if(length(reflection) < 0.1 || length(light) > sqrt(2.0)) break;
    }
    // rave_vec4 globalLightDir = {0.1,0.1,-1};
    // float global_light_participance = -dot(direction, globalLightDir);
    // // vec3 tmp = vec3(0);
    // if (global_light_participance > 0.9) {
    //     light.x += (rave_vec4){.9,.9,.6}.x * reflection.x * global_light_participance / 2.0;
    //     light.y += (rave_vec4){.9,.9,.6}.y * reflection.y * global_light_participance / 2.0;
    //     light.z += (rave_vec4){.9,.9,.6}.z * reflection.z * global_light_participance / 2.0;
    // }
    // float ndot = (dot(normal, globalLightDir));
    // float _global_light_participance = (ndot>0.1)? material.roughness * ndot * 0.5 : 0.0;
    // light += (vec3(.9,.9,.6)*3.0) * reflection * _global_light_participance / 2.0;
    // light = (rave_vec4) {1,1,1};
    return light;
}

// static int counter = 0;
static void dispatch_ray(rave_ctx* ctx, int local_x, int local_y, int local_z){
// println
    rave_vec3 _origin = rave_get_ray_pos(local_x, local_y, local_z);
    rave_vec4 origin = {
        _origin.x,
        _origin.y,
        _origin.z,
    };
    rave_vec3 _direction = rave_get_ray_dir(local_x, local_y, local_z);
    rave_vec4 direction = {
        _direction.x,
        _direction.y,
        _direction.z,
    };
// println

    // counter++;
    // random_storage = (rave_vec2){(local_x*counter)*3.23, (local_y*counter)*6.432};

    rave_vec4 light = trace_ray(ctx, origin, direction);
// println
    rave_store_light(local_x,local_y,local_z, 
    (rave_vec3){light.x, light.y, light.z}
    );
// println
}
void rave_dispatch_simd(rave_ctx* ctx, int count_x, int count_y, int count_z){
    //for now single core non-simd
    srand(666);
// println

    for(int x=0; x<count_x; x++){ 
        for(int y=0; y<count_y; y++){
            for(int z=0; z<count_z; z++){
                dispatch_ray(ctx, x, y, z);
            }
        }
    }
}

void rave_sync_simd(){

}