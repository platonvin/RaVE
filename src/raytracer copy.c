#include "raytracer.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <immintrin.h>

rave_vec2 random_storage;//TODO
float Random1D(){
    // float res = sin((dot(random_storage, vec2(12.9898, 78.233))) * 43758.5453);
    // random_storage.x = random_storage.y;
    // random_storage.y = res;
    // float _;
    // float dot = 
    //     random_storage.x * 12.9898+
    //     random_storage.y * 78.233;
    // float res = modff(sin(dot * 43758.5453123) * 43758.5453, &_);
    // random_storage.x = random_storage.y;
    // random_storage.y = res;
    float res = rand() / (float)RAND_MAX;
          res = 2.0*res - 1.0;
    return (res);
}

rave_vec3 Random3D() {
    float a,b,c;
    a = Random1D();
    b = Random1D();
    c = Random1D();
    return (rave_vec3){a,b,c};
}

rave_vec3 randomSpherePoint(rave_vec3 rand) {
    float ang1 = (rand.x + 1.0) * PI; // [-1..1) -> [0..2*PI)
    float u = rand.y; // [-1..1), cos and acos(2v-1) cancel each other out, so we arrive at [-1..1)
    float u2 = u * u;
    float sqrt1MinusU2 = sqrtf(1.0 - u2);
    float x = sqrt1MinusU2 * cosf(ang1);
    float y = sqrt1MinusU2 * sinf(ang1);
    float z = u;
    return (rave_vec3){x, y, z};
}

rave_vec3 NormalOrientedHemispherePoint(rave_vec3 rand, rave_vec3 n){
    rave_vec3 v = randomSpherePoint(rand);
    float dot = 
          v.x * n.x +
          v.y * n.y +
          v.z * n.z;
    float sign = (float) (dot > 0);
    v.x *= sign;
    v.x *= sign;
    v.x *= sign;
    return v;
}

bool initTvals(rave_vec3* tMax, rave_vec3* tDelta, rave_ivec3* blockPos, rave_vec3 rayOrigin, rave_vec3 rayDirection){
    rave_vec3 effective_origin = rayOrigin;

    rave_vec3 block_corner1;
              block_corner1.x = (floorf(effective_origin.x) - effective_origin.x)/rayDirection.x; //now corners are relative vectors
              block_corner1.y = (floorf(effective_origin.y) - effective_origin.y)/rayDirection.y; //now corners are relative vectors
              block_corner1.z = (floorf(effective_origin.z) - effective_origin.z)/rayDirection.z; //now corners are relative vectors

    rave_vec3 block_corner2;
              block_corner2.x = (floor(effective_origin.x) - effective_origin.x)/rayDirection.x  + 1.0/rayDirection.x;
              block_corner2.y = (floor(effective_origin.y) - effective_origin.y)/rayDirection.y  + 1.0/rayDirection.y;
              block_corner2.z = (floor(effective_origin.z) - effective_origin.z)/rayDirection.z  + 1.0/rayDirection.z;

    (*tMax).x = fmaxf(block_corner1.x, block_corner2.x); //1 of theese will be negative so max is just to get positive
    (*tMax).y = fmaxf(block_corner1.y, block_corner2.y);
    (*tMax).z = fmaxf(block_corner1.z, block_corner2.z);

    (*tDelta).x = 1.0 / fabsf(rayDirection.x); //how many dir vectors needeed to move 1.0 across each axys
    (*tDelta).y = 1.0 / fabsf(rayDirection.y); //how many dir vectors needeed to move 1.0 across each axys
    (*tDelta).z = 1.0 / fabsf(rayDirection.z); //how many dir vectors needeed to move 1.0 across each axys

    (*blockPos) = (rave_ivec3){effective_origin.x, effective_origin.y, effective_origin.z}; //round origin to block pos

    return true;
}

//is actually precise
bool CastRay_precise(rave_ctx* ctx, rave_vec3 rayOrigin, rave_vec3 rayDirection, 
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
    
    int current_voxel = (*ctx).rave_get_voxel(voxel_pos.x, voxel_pos.y, voxel_pos.z);

    int max_steps = (*ctx).max_steps;
    assert(max_steps != 0);

    int iterations = 0;

    // __m256
    while (true) {
        bool xLy = tMax.x <= tMax.y;
        bool xLz = tMax.x <= tMax.z;
        bool yLz = tMax.y <= tMax.z;

        //LOL no perfomance benefit currently but it was there last time i tested it TODO
        fcurrentStepDiretion.x = (int)((int)(( xLy) && ( xLz)));
        fcurrentStepDiretion.y = (int)((int)((!xLy) && ( yLz)));
        fcurrentStepDiretion.z = (int)((int)((!xLz) && (!yLz)));

        voxel_pos.x += steps.x * fcurrentStepDiretion.x;
        voxel_pos.y += steps.y * fcurrentStepDiretion.y;
        voxel_pos.z += steps.z * fcurrentStepDiretion.z;
        
        tMax.x += tDelta.x * fcurrentStepDiretion.x;
        tMax.y += tDelta.y * fcurrentStepDiretion.y;
        tMax.z += tDelta.z * fcurrentStepDiretion.z;
        
        current_voxel = (*ctx).rave_get_voxel(voxel_pos.x, voxel_pos.y, voxel_pos.z);

        // pi(voxel_pos);
        // pf(tDelta);
        // printf(voxel_pos.x)
        
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

    rave_vec3 tFinal;
              tFinal.x = tMax.x - tDelta.x;
              tFinal.y = tMax.y - tDelta.y;
              tFinal.z = tMax.z - tDelta.z;
    // block_fraction += dot(tFinal, fcurrentStepDiretion);
    block_fraction = 
        tFinal.x * (float)fcurrentStepDiretion.x+
        tFinal.y * (float)fcurrentStepDiretion.y+
        tFinal.z * (float)fcurrentStepDiretion.z;

    (*material) = (*ctx).rave_get_material(current_voxel);
    (*fraction) = block_fraction;

    return (block_hit);
}
float dot(rave_vec3 v1, rave_vec3 v2) {
    return v1.x*v2.x + 
           v1.y*v2.y + 
           v1.z*v2.z;
}
rave_vec3 reflect(rave_vec3 I, rave_vec3 N) {
    float dotProduct = dot(I, N);
    rave_vec3 result;
    result.x = I.x - 2.0f * dotProduct * N.x;
    result.y = I.y - 2.0f * dotProduct * N.y;
    result.z = I.z - 2.0f * dotProduct * N.z;
    return result;
}
rave_vec3 normalize(rave_vec3 v) {
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
rave_vec3 cross(rave_vec3 v1, rave_vec3 v2) {
    rave_vec3 result;
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}
rave_vec3 mix(rave_vec3 v1, rave_vec3 v2, float t) {
    rave_vec3 result;
    result.x = v1.x * (1.0f - t) + v2.x * t;
    result.y = v1.y * (1.0f - t) + v2.y * t;
    result.z = v1.z * (1.0f - t) + v2.z * t;
    return result;
}
float distance(rave_vec3 v1, rave_vec3 v2) {
    float dx = v1.x - v2.x;
    float dy = v1.y - v2.y;
    float dz = v1.z - v2.z;
    return sqrtf(dx * dx + dy * dy + dz * dz);
}
rave_vec3 add(rave_vec3 v1, rave_vec3 v2) {
    rave_vec3 result;
    result.x = v1.x + v2.x;
    result.y = v1.y + v2.y;
    result.z = v1.z + v2.z;
    return result;
}
rave_vec3 sub(rave_vec3 v1, rave_vec3 v2) {
    rave_vec3 result;
    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    result.z = v1.z - v2.z;
    return result;
}
rave_vec3 mul(rave_vec3 v, float scalar) {
    rave_vec3 result;
    result.x = v.x * scalar;
    result.y = v.y * scalar;
    result.z = v.z * scalar;
    return result;
}
rave_vec3 clamp(rave_vec3 v, float minVal, float maxVal) {
    rave_vec3 result;
    result.x = fminf(fmaxf(v.x, minVal), maxVal);
    result.y = fminf(fmaxf(v.y, minVal), maxVal);
    result.z = fminf(fmaxf(v.z, minVal), maxVal);
    return result;
}
void pf(rave_vec3 v) {
    printf("rave_vec3: (%f, %f, %f)\n", v.x, v.y, v.z);
}
void pi(rave_ivec3 v) {
    printf("rave_vec3: (%d, %d, %d)\n", v.x, v.y, v.z);
}
void ProcessHit(rave_vec3* origin, rave_vec3* direction, 
                float fraction, rave_vec3 normal, Material material, 
                rave_vec3* accumulated_light, rave_vec3* accumulated_reflection){

            rave_vec3 new_origin;
                      new_origin.x = (*origin).x + (fraction * (*direction).x);
                      new_origin.y = (*origin).y + (fraction * (*direction).y);
                      new_origin.z = (*origin).z + (fraction * (*direction).z);

            rave_vec3 hemisphereDistributedDirection = NormalOrientedHemispherePoint(Random3D(), normal);

            rave_vec3 new_direction = hemisphereDistributedDirection;

            rave_vec3 idealReflection = reflect((*direction), normal);
            new_direction = normalize(mix(idealReflection, new_direction, material.roughness));
            new_origin.x += normal.x * 0.001;
            new_origin.y += normal.y * 0.001;
            new_origin.z += normal.z * 0.001;

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

static rave_vec3 trace_ray(rave_ctx* ctx, rave_vec3 origin, rave_vec3 direction){
    //finds intersections and processes hits
    float fraction = 0.0;
    rave_vec3 normal = {0};
    Material material = {0};

    rave_vec3 light = {0};
    rave_vec3 reflection = {1,1,1};
// println
    for (int i = 0; (i < ctx->max_reflections); i++){
        bool hit = CastRay_precise(ctx, origin, direction, &fraction, &normal, &material);
        // if(hit) printf("yaw! %f\n", material.emmitance);
        // if(hit) printf("yaw! %f\n", fraction);
        if(!hit) break;
        
        ProcessHit(&origin, &direction, 
            fraction, normal, material, 
            &light, &reflection);
        
        // printf("%f\n", light.x);
        
        // if(length(reflection) < 0.1 || length(light) > sqrt(2.0)) break;
    }
    // rave_vec3 globalLightDir = {0.1,0.1,-1};
    // float global_light_participance = -dot(direction, globalLightDir);
    // // vec3 tmp = vec3(0);
    // if (global_light_participance > 0.9) {
    //     light.x += (rave_vec3){.9,.9,.6}.x * reflection.x * global_light_participance / 2.0;
    //     light.y += (rave_vec3){.9,.9,.6}.y * reflection.y * global_light_participance / 2.0;
    //     light.z += (rave_vec3){.9,.9,.6}.z * reflection.z * global_light_participance / 2.0;
    // }
    // float ndot = (dot(normal, globalLightDir));
    // float _global_light_participance = (ndot>0.1)? material.roughness * ndot * 0.5 : 0.0;
    // light += (vec3(.9,.9,.6)*3.0) * reflection * _global_light_participance / 2.0;
    // light = (rave_vec3) {1,1,1};
    return light;
}

// static int counter = 0;
static void dispatch_ray(rave_ctx* ctx, int local_x, int local_y, int local_z){
// println
    rave_vec3 origin = ctx->rave_get_ray_pos(local_x, local_y, local_z);
// println
    rave_vec3 direction = ctx->rave_get_ray_dir(local_x, local_y, local_z);
// println

    // counter++;
    // random_storage = (rave_vec2){(local_x*counter)*3.23, (local_y*counter)*6.432};

    rave_vec3 light = trace_ray(ctx, origin, direction);
// println
    ctx->rave_store_light(local_x,local_y,local_z, light);
// println
}
void rave_dispatch(rave_ctx* ctx, int count_x, int count_y, int count_z){
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

void rave_sync(){

}