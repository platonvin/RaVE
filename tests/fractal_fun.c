#define RAVE_CUSTOM_VOXEL_TYPE int

#include "../src/raytracer.c"

#include <assert.h>
#include <stdint.h>

float dot4(rave_vec4 q1, rave_vec4 q2){
    return 
        q1.x*q2.x+
        q1.y*q2.y+
        q1.z*q2.z+
        q1.w*q2.w;
}
rave_vec4 quat_mul(rave_vec4 q1, rave_vec4 q2){
    rave_vec4 r;
    rave_vec3 q1_yzw = (rave_vec3){q1.y, q1.z, q1.w};
    rave_vec3 q2_yzw = (rave_vec3){q2.y, q2.z, q2.w};
    
    r.x = q1.x*q2.x - rave_dot3(q1_yzw, q2_yzw);


    // r.yzw = q1.x*q2.yzw + q2.x*q1.yzw + cross3( q1.yzw, q2.yzw );
    r.y = rave_mul3(q2_yzw, q1.x).x + rave_mul3(q1_yzw, q2.x).x + rave_cross3(q1_yzw, q2_yzw).x;
    r.z = rave_mul3(q2_yzw, q1.x).y + rave_mul3(q1_yzw, q2.x).y + rave_cross3(q1_yzw, q2_yzw).y;
    r.w = rave_mul3(q2_yzw, q1.x).z + rave_mul3(q1_yzw, q2.x).z + rave_cross3(q1_yzw, q2_yzw).z;
    
    return r;
}
rave_vec4 quat_square(rave_vec4 q){
    rave_vec4 r;
    rave_vec3 q_yzw = (rave_vec3){q.y, q.z, q.w};

    r.x = q.x*q.x - rave_dot3(q_yzw, q_yzw);

    // r.yzw = 2.0*q.x*q.yzw;
    r.y = 2.0*q.x*q.y;
    r.z = 2.0*q.x*q.z;
    r.w = 2.0*q.x*q.w;
    return r;
}
float quat_length_squared(rave_vec4 q){
    return dot4(q,q);
}
rave_vec4 quat_cube(rave_vec4 q){
    rave_vec4 q2;
        q2.x = q.x*q.x;
        q2.y = q.y*q.y;
        q2.z = q.z*q.z;
        q2.w = q.w*q.w;
    
    return (rave_vec4){
        q.x * (q2.x - 3.0*q2.y - 3.0*q2.z - 3.0*q2.w), 
        q.y * (3.0*q2.x - q2.y - q2.z - q2.w),
        q.z * (3.0*q2.x - q2.y - q2.z - q2.w),
        q.w * (3.0*q2.x - q2.y - q2.z - q2.w)};
}

const float   ESCAPE_THRESHOLD = 20.0;
const float CONVERGE_THRESHOLD = 0.3;
const int   MAX_ITERATIONS = 12; //actually, very few steps needed

float check_julia(rave_vec3 pos, rave_vec4 C, int* converged_i){
    rave_vec4 z = (rave_vec4){pos.x, pos.y, pos.z, 0.0};
	float m2 = 0.0;
    bool did_converge = false;

    for(int i=0; i < MAX_ITERATIONS; i++) {
        //for ^3
		// z = quat_cube( z ) + C;
        // printf("%f:%f:%f:%f\n", z.x,z.y,z.z,z.w);
        rave_vec4 z3 = quat_cube(z);
		z.x = z3.x + C.x;
		z.y = z3.y + C.y;
		z.z = z3.z + C.z;
		z.w = z3.w + C.w;
        // printf("%f:%f:%f:%f\n", z.x,z.y,z.z,z.w);

        //for ^2
		// z = quat_square( z ) + C;
        
        m2 = quat_length_squared(z);
        // printf("%d:%f\n", i, m2);
        if(m2 > ESCAPE_THRESHOLD) {
            (*converged_i) = i;
            break;
        }
        if((m2 < CONVERGE_THRESHOLD) && (!did_converge)){
            (*converged_i) = i;
            break;
        }
	}

	return m2;        
}