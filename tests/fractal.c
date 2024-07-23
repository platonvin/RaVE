//julia set

#include "../dependencies/tgafunc.h"
// #include "../src/raytracer_simd.h"

// #include <assert.h>
// #include <stdint.h>
#include <math.h>
#include <stdio.h>
// #include <stdlib.h>

#include "fractal_fun.c"

void delete_if_exists(const char *filename) {
    assert(filename != NULL);
    FILE *file = fopen(filename, "r");
    assert((file && !fclose(file) && !remove(filename)) || !file);
}
voxel get_bounding_box_voxel(int x, int y, int z){
    if((x==-70) || (y==-70) || (z==-70) ||
       (x==+70) || (y==+70) || (z==+70)
    ) return 1;
    else 
    return 0;
}
voxel get_light_voxel(int x, int y, int z){
    if(z==69) return 2;
    else return 0;
}
voxel rave_get_voxel(int x, int y, int z){
    // C_x = -0.65 + -0.65*0.15 * Math.sin((currentTime + 23.0)/2.12);
    // C_y = -0.3  + -0.3 *0.15 * Math.sin((currentTime + 23.0)/3.523);
    // C_z = +0.6  + +0.6 *0.15 * Math.sin((currentTime + 23.0)/5.634);
    // C_w = -0.2  + -0.2 *0.15 * Math.sin((currentTime + 23.0)/7.6345);
    
    rave_vec4 Constant = {
        -0.45, 
        +0.25, 
        -0.60, 
        +0.30,
        };
    rave_vec3 pos = {
        (((float)x) + 0.5f) / 10.0,
        (((float)y) + 0.5f) / 10.0,
        (((float)z) + 0.5f) / 10.0};

    // printf("pos %f %f %f\n", pos.x, pos.y, pos.z);
    // if (z < 0) {
        float julia = check_julia(pos, Constant);
        if (julia < ESCAPE_THRESHOLD) {
            // printf("js: %f\n", julia);
            float division = (0.5+sin(julia * 5624536.0)*0.5)*10.0;
            return (int)division + 1;
        }
    // }
    return get_bounding_box_voxel(x,y,z) + get_light_voxel(x,y,z);
}

static int width = 512, height = 512;
uint32_t max_reflections = 4;
uint32_t max_steps = 256;
const int sample_count = 4;

static rave_vec3 camera_ray_dir;
static rave_vec3 camera_ray_dir_plane;
static rave_vec3 horizline;
static rave_vec3 vertiline;

rave_vec3 rave_get_ray_pos(int x, int y, int z){
    rave_vec3 pos = {1.0*25.0, 0.4*25.0, 0.5*25.0};

    float uv_x = (x*2 - width ) / (float)width ;
    float uv_y = (y*2 - height) / (float)height;

    float scale_x = uv_x * 35.0;
    float scale_y = uv_y * 35.0;
    
    pos = add3(pos, mul3(horizline, scale_x));
    pos = add3(pos, mul3(vertiline, scale_y));

    return pos;
}
rave_vec3 rave_get_ray_dir(int x, int y, int z){
    return camera_ray_dir;
}


static uint8_t *data;
static tga_info *info;
void rave_store_light(int x, int y, int z, rave_vec3 light){
    uint8_t *pixel;
    rave_vec3 final_light = clamp3(light, 0, 1);
    pixel = tga_get_pixel(data, info, x, y);
    pixel[0] += (uint8_t)(final_light.x * (255.0 / ((float)sample_count)));
    pixel[1] += (uint8_t)(final_light.y * (255.0 / ((float)sample_count)));
    pixel[2] += (uint8_t)(final_light.z * (255.0 / ((float)sample_count)));
}
Material rave_get_material(voxel voxel){
    Material mat = {0};
    mat.color = (rave_vec3){0.65,.725,.8};
    mat.emmitance = 0.01;
    mat.roughness = 0.95;
    if (voxel == 2) {
        mat.color = (rave_vec3){0.65,.725,.8};
        mat.emmitance = 0.9;
    } else {
        mat.color.y = sin(((float) voxel)*5152.0)*0.5 + 0.5;
        mat.emmitance = ((float) voxel) / 20.0;
    }
    // if (voxel == 3) {
    //     mat.color = (rave_vec3){0.1,.9,.1};
    //     mat.emmitance = 0.3;
    // }
    return mat;
}
int main(int argc, char *argv[]) {
    const char* out_name = "image.tga";

    camera_ray_dir = normalize3((rave_vec3){-1.0, -0.4, -0.5});
    camera_ray_dir_plane = normalize3((rave_vec3){camera_ray_dir.x, camera_ray_dir.y, 0});
    horizline = normalize3(cross3(camera_ray_dir_plane, (rave_vec3){0,0,1}));
    // horizline = normalize((rave_vec3){1.0f, -1.0f, 0.0f});
    vertiline = normalize3(cross3(camera_ray_dir, horizline));
    // globalLightDir = normalize(mul((rave_vec3){12.0f, 21.0f, 7.0f}, -1.0f));

    delete_if_exists(out_name);
    enum tga_error error_code;
    error_code = tga_create(&data, &info, width, height, TGA_PIXEL_RGB24);
    assert(error_code == TGA_NO_ERROR);

    rave_ctx rave_state = {};
    rave_init();

// println
    rave_dispatch_simd(&rave_state, width, height, sample_count);    

    // Saves the image as a TGA file.
    error_code = tga_save_from_info(data, info, out_name);
    assert(error_code == TGA_NO_ERROR);

    tga_free_data(data);
    tga_free_info(info);

    return 0;
}