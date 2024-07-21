#include "../dependencies/tgafunc.h"
#include "../src/raytracer.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// #include <immintrin.h>
// #include <unistd.h>
// #include "../src/raytracer.h"
void delete_if_exists(const char *filename) {
    assert(filename != NULL);
    FILE *file = fopen(filename, "r");
    assert((file && !fclose(file) && !remove(filename)) || !file);
}
voxel get_box_voxel(int x, int y, int z){
    if((x==-20) || (y==-20) || (z==-20) ||
       (x==+20) || (y==+20) || (z==+20)
    ) return 1;
    else return 0;
}
voxel get_sphere_voxel(int x, int y, int z){
    rave_vec3 sphere_center = (rave_vec3) {0,0,0};
    rave_vec3 vox_center = add3((rave_vec3){x,y,z}, (rave_vec3){0.5,0.5,0.5});
    float d = distance3(vox_center, sphere_center);
    // return 0;
    if(d < 2.0){
        return 3;
    } else return 0;
}
voxel get_light_voxel(int x, int y, int z){
    if(z==19) return 2;
    else return 0;
}
voxel rave_get_voxel(int x, int y, int z){
    // rave_vec3 vox_center = (rave_vec3){x,y,z};
    voxel box_voxel = get_box_voxel(x,y,z);
    voxel sphere_voxel = get_sphere_voxel(x,y,z);
    voxel light_voxel = get_light_voxel(x,y,z);

    return box_voxel + sphere_voxel + light_voxel;
}

static int width = 512, height = 512;

static rave_vec3 cameraRayDir;
static rave_vec3 cameraRayDirPlane;
static rave_vec3 horizline;
static rave_vec3 vertiline;
static rave_vec3 cameraPos;
static rave_vec3 globalLightDir;
rave_vec3 rave_get_ray_pos(int x, int y, int z){
    rave_vec3 pos = {9.002,7.003,6.504};

    float uv_x = (x*2 - width ) / (float)width ;
    float uv_y = (y*2 - height) / (float)height;

    float scale_x = uv_x * 5.0;
    float scale_y = uv_y * 5.0;
    
    pos = add3(pos, mul3(horizline, scale_x));
    pos = add3(pos, mul3(vertiline, scale_y));

    return pos;
}
rave_vec3 rave_get_ray_dir(int x, int y, int z){
    return cameraRayDir;
}


static uint8_t *data;
static tga_info *info;
void rave_store_light(int x, int y, int z, rave_vec3 light){
    uint8_t *pixel;
    rave_vec3 final_light = clamp3(light, 0, 1);
    pixel = tga_get_pixel(data, info, x, y);
    pixel[0] += (uint8_t)(final_light.x * (255.0 / 5.0));
    pixel[1] += (uint8_t)(final_light.y * (255.0 / 5.0));
    pixel[2] += (uint8_t)(final_light.z * (255.0 / 5.0));
}
Material rave_get_material(voxel voxel){
    Material mat = {0};
    mat.color = (rave_vec3){0.3,.4,.9};
    mat.emmitance = 0.01;
    mat.roughness = 0.4;
    if (voxel == 2) {
        mat.color = (rave_vec3){0.1,.2,.9};
        mat.emmitance = 0.9;
    }
    // if (voxel == 3) {
    //     mat.color = (rave_vec3){0.1,.9,.1};
    //     mat.emmitance = 0.3;
    // }
    return mat;
}
int main(int argc, char *argv[]) {
    const char* out_name = "image.tga";

    cameraRayDir = normalize3((rave_vec3){-1.0, -0.8, -0.6});
    cameraRayDirPlane = normalize3((rave_vec3){cameraRayDir.x, cameraRayDir.y, 0});
    horizline = normalize3(cross3(cameraRayDirPlane, (rave_vec3){0,0,1}));
    // horizline = normalize((rave_vec3){1.0f, -1.0f, 0.0f});
    vertiline = normalize3(cross3(cameraRayDir, horizline));
    // globalLightDir = normalize(mul((rave_vec3){12.0f, 21.0f, 7.0f}, -1.0f));

    delete_if_exists(out_name);
    enum tga_error error_code;
    error_code = tga_create(&data, &info, width, height, TGA_PIXEL_RGB24);
    assert(error_code == TGA_NO_ERROR);

    rave_ctx rave_state = {};
    rave_state.max_steps = 512;
    rave_state.max_reflections = 5;
    // rave_state.rave_get_material
    rave_state.rave_get_ray_dir  = &rave_get_ray_dir;
    rave_state.rave_get_material = &rave_get_material;
    rave_state.rave_get_ray_pos  = &rave_get_ray_pos;
    rave_state.rave_get_voxel    = &rave_get_voxel;
    rave_state.rave_store_light  = &rave_store_light;
println
    rave_dispatch(&rave_state, width, height, 5);    

    // Saves the image as a TGA file.
    error_code = tga_save_from_info(data, info, out_name);
    assert(error_code == TGA_NO_ERROR);

    tga_free_data(data);
    tga_free_info(info);

    return 0;
}