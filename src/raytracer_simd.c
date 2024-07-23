#include "raytracer_simd.h"
#include <pthread.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <immintrin.h>
// #include <emmintrin.h>
#include "count_cores.c"

static __m128i _my_and_epi32(__m128i __a, __m128i __b)
{
  return (__m128i)((__v4su)__a & (__v4su)__b);
}

// -1..+1
// rave_vec4 Random4D() {
//     // float a,b,c;
//     // a = Random1D();
//     // b = Random1D();
//     // c = Random1D();
//     __m128i rnd = rand_simd();
    
    
//     //Float Bits as Int
//     #define fbi(x) (*(int*)&x)

    
//     __m128i frame = (
//         _my_and_epi32(rnd, 
//         _mm_set1_epi32(0b111111111111111)));
//     __m128 frnd = _mm_sub_ps(
//         _mm_div_ps(
//             _mm_cvtepi32_ps((frame)), 
//             _mm_set1_ps(((float)0b111111111111111))/2.0),
//         _mm_set1_ps(1.0));

//     float* ptr = &frnd;
//     // sum = _mm_add_ps(sum, frnd);
//     // pf(frnd);
//     return (rave_vec4){ptr[0], ptr[1], ptr[2], ptr[3]};
//     // return (rave_vec4){a,b,c,1.0};
// }
float hsum_ps_sse3(__m128 v) {
    __m128 shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
    __m128 sums = _mm_add_ps(v, shuf);
    shuf        = _mm_movehl_ps(shuf, sums); // high half -> low half
    sums        = _mm_add_ss(sums, shuf);
    return        _mm_cvtss_f32(sums);
}

rave_vec3 random_sphere_point_simd() {
    float r1 = ((float)rand()) / ((float) RAND_MAX)*2.0 - 1.0;
    float r2 = ((float)rand()) / ((float) RAND_MAX)*2.0 - 1.0;
    
    float ang1 = (r1 + 1.0) * PI; // [-1..1) -> [0..2*PI)
    float u = r2; // [-1..1), cos and acos(2v-1) cancel each other out, so we arrive at [-1..1)
    float u2 = u * u;
    float sqrt1MinusU2 = sqrtf(1.0 - u2);
    float x = sqrt1MinusU2 * cosf(ang1);
    float y = sqrt1MinusU2 * sinf(ang1);
    float z = u;
    return (rave_vec3){x, y, z};
}

__m128 normal_oriented_hemisphere_point_simd(__m128 normal){
    rave_vec3 v = random_sphere_point_simd();
    __m128 r = _mm_setr_ps(v.x, v.y, v.z, 0);
    float dot = 
          v.x * normal[0] +
          v.y * normal[1] +
          v.z * normal[2];
    float sign = (dot>0) ? 1.0 : -1.0;
    r = _mm_mul_ps(r, _mm_set1_ps(sign));

    return r;
}

rave_vec3 normal_oriented_hemisphere_point(rave_vec3 n){
    rave_vec3 v = random_sphere_point_simd();
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
inline bool initTvals_simd(__m128* restrict tMax, __m128* restrict tDelta, __m128i* restrict blockPos, __m128 rayOrigin, __m128 rayDirection){


    __m128 effective_origin = rayOrigin;
    __m128 ray_direction = rayDirection;

    __m128 block_corner1 = 
        _mm_div_ps(
            _mm_sub_ps(
                _mm_floor_ps(effective_origin),
                effective_origin
            ),
            ray_direction
        );
    __m128 block_corner2 = 
        _mm_add_ps(
            block_corner1,
            _mm_div_ps(
                _mm_set1_ps(1.0),
                ray_direction
            )
        );

    (*tMax) = _mm_max_ps(block_corner1, block_corner2);

    (*tDelta) = _mm_div_ps(
        _mm_set1_ps(1.0),
        _mm_abs_ps(ray_direction));

    //_MM_FROUND_FLOOR  / _MM_FROUND_CEIL should be used instead of _MM_FROUND_TRUNC
    (*blockPos) = _mm_cvtps_epi32(_mm_round_ps((effective_origin), _MM_FROUND_NO_EXC | _MM_FROUND_TO_NEG_INF));

    return true;
}
//is actually precise
int CastRay_precise_simd(__m128 rayOrigin, __m128 rayDirection, 
        float* fraction, __m128* normal, Material* material){
    bool block_hit = false;

    __m128i steps = {};

    steps = _mm_castps_si128(_mm_cmpgt_ps(
        // _mm_set_ps(rayDirection.x, rayDirection.y, rayDirection.z, rayDirection.w), 
        rayDirection, 
        _mm_set1_ps(0)));

    steps = _my_and_epi32(steps, _mm_set1_epi32(1));
    steps = _mm_mullo_epi32(steps, _mm_set1_epi32(2));
    steps = _mm_sub_epi32(steps, _mm_set1_epi32(1));

    // rave_vec4 fsteps = {steps.x, steps.y, steps.z};
    // __m128 fsteps = _mm_ps
    // rave_vec4 fsteps;
    // fsteps.vec = _mm_cvtepi32_ps(steps.vec);

    __m128 tMax = {0};
    __m128 tDelta = {0};
    __m128i voxel_pos = {0};
    float block_fraction = 0.0;

    // __m128 fcurrentStepDiretion = {0};

    initTvals_simd(&tMax, &tDelta, &voxel_pos, rayOrigin, rayDirection); //does not intersect with scene
    // pf(tMax);
    // pf(tDelta);
    // pi(voxel_pos);

    // initTvals(&tMax, &tDelta, &voxel_pos, (rave_vec3){rayOrigin.x, rayOrigin.y, rayOrigin.z}, (rave_vec3){rayDirection.x, rayDirection.y, rayDirection.z}); //does not intersect with scene
    
    // _mm_storeu_ps(&_vp, _mm_castsi128_ps(voxel_pos));
    // printf("_vp %d %d %d %d\n", _vp[0], _vp[1], _vp[2], _vp[3]);
    int* _vp = (int*) &voxel_pos;
    int current_voxel = rave_get_voxel(_vp[0], _vp[1], _vp[2]);

    int iterations = 0;

    // __m256
    __m128i step_dir = {};
    while (true) {

        __m128 tMax_xyzw = tMax;
        __m128 tMax_yzxw = _mm_shuffle_ps(tMax_xyzw, tMax_xyzw, _MM_SHUFFLE(3,0,2,1));

        __m128 compared = _mm_cmple_ps(tMax_xyzw, tMax_yzxw);
        __m128i Less_XY_YZ_ZX = _mm_castps_si128(compared);
        __m128i Less_ZX_XY_YZ = _mm_castps_si128(
            _mm_shuffle_ps(
                _mm_castsi128_ps(Less_XY_YZ_ZX), 
                _mm_castsi128_ps(Less_XY_YZ_ZX), 
                _MM_SHUFFLE(3,1,0,2)
            )
        );
        __m128i Greater_ZX_XY_YZ = _mm_castps_si128(_mm_xor_ps(_mm_castsi128_ps(Less_ZX_XY_YZ), _mm_castsi128_ps(_mm_set1_epi32(0xffffffff))));
        step_dir = _my_and_epi32(_my_and_epi32(Less_XY_YZ_ZX, Greater_ZX_XY_YZ), _mm_set_epi32(1,1,1,1));

        // pi(voxel_pos)
        // pi(steps)
        // pi(step_dir)
        
        voxel_pos = _mm_add_epi32(
            voxel_pos, 
            _mm_mullo_epi32(steps, step_dir)
        );
        // pi(voxel_pos)

        __m128 fstep_dir = _mm_cvtepi32_ps(step_dir);
        tMax = _mm_add_ps(
            tMax, 
            _mm_mul_ps(tDelta, fstep_dir)
        );

        int* _vp = (int*) &voxel_pos;
        current_voxel = rave_get_voxel(_vp[0], _vp[1], _vp[2]);

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

    __m128 _normal = _mm_cvtepi32_ps(_mm_mullo_epi32(_mm_mullo_epi32(steps, step_dir), _mm_set1_epi32(-1)));
    // float* _nm = (float*) &_normal;
        // _mm_storeu_ps(_nm, _normal);
        (*normal) = _normal;

    __m128 tFinal = _mm_sub_ps(tMax, tDelta);
    __m128 _tFinal_X_step_dir = _mm_mul_ps(tFinal, _mm_cvtepi32_ps(step_dir));

    float* _bf = (float*) &_tFinal_X_step_dir;
    block_fraction = _bf[0] + _bf[1] + _bf[2];

    // printf("simd_voxel %d \n", current_voxel);
    (*material) = rave_get_material(current_voxel);
    (*fraction) = block_fraction;

    return (block_hit);
    // return (block_hit);
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

// float dot(rave_vec4 v1, rave_vec4 v2) {
//     return v1.x*v2.x + 
//            v1.y*v2.y + 
//            v1.z*v2.z;
// }
__m128 reflect_simd(__m128 I, __m128 N) {
    float dotProduct;
    dotProduct = 
        I[3]*N[3]+
        I[2]*N[2]+
        I[1]*N[1];
    
    __m128 result = _mm_sub_ps(
        I,
        _mm_mul_ps(N, _mm_set1_ps(2.0*dotProduct))
    );
    // result.x = I.x - 2.0f * dotProduct * N.x;
    // result.y = I.y - 2.0f * dotProduct * N.y;
    // result.z = I.z - 2.0f * dotProduct * N.z;
    return result;
}
__m128 normalize(__m128 v) {
    float length = sqrtf(v[3]*v[3]+ 
                         v[2]*v[2]+ 
                         v[1]*v[1]);
                         
    __m128 result = {};

    assert(length != 0.0f);
    assert(length > 0.0f);
    if (length > 0.0f) {
        // result.x = v.x / length;
        // result.y = v.y / length;
        // result.z = v.z / length;
        result = _mm_div_ps(v, _mm_set1_ps(length));
    } else {
        // abort();
        // result = (rave_vec4){};
    }
    return result;
}
// rave_vec4 cross(rave_vec4 v1, rave_vec4 v2) {
//     rave_vec4 result;
//     result.x = v1.y * v2.z - v1.z * v2.y;
//     result.y = v1.z * v2.x - v1.x * v2.z;
//     result.z = v1.x * v2.y - v1.y * v2.x;
//     return result;
// }
__m128 mix(__m128 v1, __m128 v2, float t) {
    __m128 result;
    assert(t >= 0.0);
    assert(t <= 1.0);

    // result.x = v1.x * (1.0f - t) + v2.x * t;
    // result.y = v1.y * (1.0f - t) + v2.y * t;
    // result.z = v1.z * (1.0f - t) + v2.z * t;
    result = _mm_add_ps(
        _mm_mul_ps(v1, _mm_set1_ps(1.0f-t)),
        _mm_mul_ps(v2, _mm_set1_ps(t))
    );
    // result.w = v1.w * (1.0f - t) + v2.w * t;
    return result;
}
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
void ProcessHit_simd(__m128* origin, __m128* direction, 
                float fraction, __m128 normal, Material material, 
                __m128* accumulated_light, __m128* accumulated_reflection){

            // rave_vec4 new_origin;
            //           new_origin.x = (*origin).x + (fraction * (*direction).x);
            //           new_origin.y = (*origin).y + (fraction * (*direction).y);
            //           new_origin.z = (*origin).z + (fraction * (*direction).z);

            __m128 new_origin = _mm_add_ps(_mm_mul_ps(*direction, _mm_set1_ps(fraction)), *origin);
// println
            // normal = normalize(normal);

            // rave_vec4 rnd = Random4D();
            
            __m128 hemisphereDistributedDirection = normal_oriented_hemisphere_point_simd(normal);
// println 

                // printf("rnd->x %f\n", rnd.x);
                // printf("rnd->y %f\n", rnd.y);
                // printf("rnd->z %f\n", rnd.z);

                // printf("normal->x %f\n", normal.x);
                // printf("normal->y %f\n", normal.y);
                // printf("normal->z %f\n", normal.z);

                // printf("hemisphereDistributedDirection->x %f\n", hemisphereDistributedDirection.x);
                // printf("hemisphereDistributedDirection->y %f\n", hemisphereDistributedDirection.y);
                // printf("hemisphereDistributedDirection->z %f\n", hemisphereDistributedDirection.z);

            // rave_vec4 new_direction = normalize((rave_vec4){hemisphereDistributedDirection[3], hemisphereDistributedDirection[2], hemisphereDistributedDirection[1]});
            // rave_vec4 new_direction = (hemisphereDistributedDirection);

            // bool ass = fabsf(normal.x + normal.y + normal.z) == 1.0;
            // if(!ass){
            //     printf("normal->x %f\n", normal.x);
            //     printf("normal->y %f\n", normal.y);
            //     printf("normal->z %f\n", normal.z);
            //     printf("\n");
            //     abort();
            // }
            // assert();

            // rave_vec4 idealReflection = (reflect((*direction), normal));
// println
            __m128 idealReflection = normalize(reflect_simd((*direction), normal));

// println
            __m128 new_direction = normalize(mix(idealReflection, hemisphereDistributedDirection, material.roughness));
            // new_origin.x += normal.x * 0.01;
            // new_origin.y += normal.y * 0.01;
            // new_origin.z += normal.z * 0.01;
            new_origin = _mm_add_ps(new_origin, _mm_mul_ps(normal, _mm_set1_ps(0.001)));

            // (*accumulated_reflection).x *= material.color.x;
            // (*accumulated_reflection).y *= material.color.y;
            // (*accumulated_reflection).z *= material.color.z;
            __m128 color = _mm_setr_ps(material.color.x, material.color.y, material.color.z, 0);
            (*accumulated_reflection) = _mm_mul_ps((*accumulated_reflection), color);
            
            // (*accumulated_light).x += (*accumulated_reflection).x * material.emmitance;
            // (*accumulated_light).y += (*accumulated_reflection).y * material.emmitance;
            // (*accumulated_light).z += (*accumulated_reflection).z * material.emmitance;
            (*accumulated_light) = _mm_add_ps((*accumulated_light), _mm_mul_ps((*accumulated_reflection), _mm_set1_ps(material.emmitance)));
            // accumulated_light += vec3(0.8) * 0.1 * accumulated_reflection;

            // direction = reflect(direction,normal);
            (*direction) = new_direction;
            (*origin) = new_origin;

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

    (*blockPos) = (rave_ivec3){(effective_origin.x), (effective_origin.y), (effective_origin.z)}; //round origin to block pos

    return true;
}

int CastRay_precise(rave_vec3 rayOrigin, rave_vec3 rayDirection, 
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

    return (block_hit);
}
void ProcessHit(rave_vec3* origin, rave_vec3* direction, 
                float fraction, rave_vec3 normal, Material material, 
                rave_vec3* accumulated_light, rave_vec3* accumulated_reflection){

            rave_vec3 new_origin;
                      new_origin.x = (*origin).x + ((fraction) * (*direction).x);
                      new_origin.y = (*origin).y + ((fraction) * (*direction).y);
                      new_origin.z = (*origin).z + ((fraction) * (*direction).z);
// println
            normal = normalize3(normal);

            
            rave_vec3 hemisphereDistributedDirection = normal_oriented_hemisphere_point(normal);
// println 

            rave_vec3 new_direction = normalize3(hemisphereDistributedDirection);
            // rave_vec4 new_direction = (hemisphereDistributedDirection);

            // bool ass = fabsf(normal.x + normal.y + normal.z) == 1.0;
            // if(!ass){
            //     printf("normal->x %f\n", normal.x);
            //     printf("normal->y %f\n", normal.y);
            //     printf("normal->z %f\n", normal.z);
            //     printf("\n");
            //     abort();
            // }
            // assert();

            // rave_vec4 idealReflection = (reflect((*direction), normal));
// println
            rave_vec3 idealReflection = normalize3(reflect3((*direction), normal));

// println
            new_direction = normalize3(mix3(idealReflection, new_direction, material.roughness));
            new_origin.x += normal.x * 0.001f;
            new_origin.y += normal.y * 0.001f;
            new_origin.z += normal.z * 0.001f;

            (*accumulated_reflection).x *= material.color.x;
            (*accumulated_reflection).y *= material.color.y;
            (*accumulated_reflection).z *= material.color.z;
            
            (*accumulated_light).x += (*accumulated_reflection).x * material.emmitance;
            (*accumulated_light).y += (*accumulated_reflection).y * material.emmitance;
            (*accumulated_light).z += (*accumulated_reflection).z * material.emmitance;

            (*direction) = new_direction;
            (*origin) = new_origin;

}

static __m128 trace_ray_simd(__m128 origin, __m128 direction){
    //finds intersections and processes hits
    float fraction = 0.0;
    __m128 normal = {0};
    Material material = {0};

    __m128 light = {0};
    __m128 reflection = {1,1,1,1};

    for (int i = 0; (i < max_reflections); i++){
        
        int hit_s=1;
        hit_s = CastRay_precise_simd(origin, direction, &fraction, &normal, &material);
        if(!hit_s) break;

        ProcessHit_simd(&origin, &direction, 
            fraction, normal, material, 
            &light, &reflection);
        // light = _mm_set1_ps((float)hit_s / 256.0);
        
        if(hsum_ps_sse3(reflection) < 0.05 || hsum_ps_sse3(light) > (2.5)) break;
    }
    return light;
}

static rave_vec3 trace_ray(rave_vec3 origin, rave_vec3 direction){
    //finds intersections and processes hits
    float fraction = 0.0;
    rave_vec3 normal = {0};
    Material material = {0};

    rave_vec3 light = {0};
    rave_vec3 reflection = {1,1,1};
// println
    for (int i = 0; (i < max_reflections); i++){
        
        int hit=1;
        hit   = CastRay_precise(origin, direction, &fraction, &normal, &material);
        // printf("fraction= %f\n", fraction);
        // printf("material= %f %f %f _ %f %f\n", material.color.x, material.color.y, material.color.z, material.emmitance, material.roughness);
        // printf("normal= %f %f %f %f\n", normal.x, normal.y, normal.z, normal.w);
        // printf("hit= %d\n", hit);
        // int hit_s=1;
        // hit_s = CastRay_precise_simd(ctx, origin, direction, &fraction, &normal, &material);
        // printf("SIMD fraction= %f\n", fraction);
        // printf("SIMD material= %f %f %f _ %f %f\n", material.color.x, material.color.y, material.color.z, material.emmitance, material.roughness);
        // printf("SIMD normal= %f %f %f %f\n", normal.x, normal.y, normal.z, normal.w);
        // printf("SIMD hit= %d\n", hit_s);
        // printf("\n");
        if(!hit) break;

        // if(!hit_s) break;

        ProcessHit(&origin, &direction, 
            fraction, normal, material, 
            &light, &reflection);
        // if(length(reflection) < 0.1 || length(light) > sqrt(2.0)) break;
    }
    return light;
}

// static int counter = 0;
void rave_init(){
    srand(69);
    // next_simd = _mm_set_epi32(42,888888888,69,420);
}
typedef struct dispatch_args_t{
    int initial;
    int total;
    int step;
    int total_x;
    int total_y;
    int total_z;
} dispatch_args_t;

static void dispatch_ray(int local_x, int local_y, int local_z){
// println 
    rave_vec3 _origin = rave_get_ray_pos(local_x, local_y, local_z);
    __m128 origin = _mm_setr_ps(
        _origin.x,
        _origin.y,
        _origin.z,
        0
    );
    rave_vec3 _direction = rave_get_ray_dir(local_x, local_y, local_z);
    __m128 direction = _mm_setr_ps(
        _direction.x,
        _direction.y,
        _direction.z,
        0
    );


    __m128 light = trace_ray_simd(origin, direction);
    rave_store_light(local_x, local_y, local_z, 
    (rave_vec3){light[0], light[1], light[2]}
    );
    // rave_vec3 light = trace_ray(_origin, _direction);
    // rave_store_light(local_x, local_y, local_z, 
    // (rave_vec3){light.x, light.y, light.z}
    // );
}

static void* dispatch_ray_wave(void* _args){
    dispatch_args_t* args = _args;
    srand((*args).initial*42);
    // single core simd
    // for(int x=0; x<count_x; x++){ 
    //     for(int y=0; y<count_y; y++){
    //         for(int z=0; z<count_z; z++){
    //             dispatch_ray(ctx, x, y, z);
    //         }
    //     }
    // }
    // printf()
    // printf("(*args).initial = %d\n", (*args).initial);
    // printf("(*args).total = %d\n", (*args).total);
    // printf("(*args).step = %d\n", (*args).step);
    // (*args).total_y /= 4;
    // (*args).total_y *= 1+(*args).initial;
    for(int i=(*args).initial; i<(*args).total; i+=(*args).step){
        int x =  i                     % (*args).total_x ;
        int y = (i /  (*args).total_x) % (*args).total_y ;
        int z =  i / ((*args).total_x  * (*args).total_y);

        // printf("%d %d %d\n", x, y, z);
        dispatch_ray(x, y, z);
    }
    
    return NULL;
}

void rave_dispatch_simd(rave_ctx* ctx, int count_x, int count_y, int count_z){
    int thread_count = get_core_count();
    // thread_count = 1;
    
    assert(thread_count != 0);
    pthread_t* threads = malloc(thread_count * sizeof(pthread_t));
    
    for (int i=0; i<thread_count; i++) {
        dispatch_args_t* args = malloc(sizeof(dispatch_args_t)); // Allocate memory for arguments

        (*args).initial = i;
        (*args).step = thread_count;
        (*args).total = (count_x*count_y*count_z);
        (*args).total_x = count_x;
        (*args).total_y = count_y;
        (*args).total_z = count_z;
        
        int pthread_creation = pthread_create(&threads[i], NULL, dispatch_ray_wave, (void*)args);

        assert(pthread_creation == 0);
    }
    for (int i=0; i<thread_count; i++) {
        int pthread_joining = pthread_join(threads[i], NULL);
        assert(pthread_joining == 0);
    }

    free(threads);
}

void rave_sync_simd(){

}