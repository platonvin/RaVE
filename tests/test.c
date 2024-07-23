#include "../dependencies/tgafunc.h"
#include "../src/raytracer_simd.h"

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

enum voxel_e {
    EMPTY,
    RED_WALL,
    GREEN_WALL,
    WHITE_WALL,
    GLOSSY_WALL,
    LIGHT,
    SPHERE,
    BOX,
};

voxel get_bounding_box_voxel(int x, int y, int z){
    if(x==-30) {
        return RED_WALL;
    } else 
    if (x==+30) {
        return GREEN_WALL;
    } else 
    if(z==+30){
        return GLOSSY_WALL;
    } else
    if ((x<=-30) || (y<=-30) || (z<=-30) || 
        (x>=+30) || (y>=+30) || (z>=+30)) {
        return WHITE_WALL;
    }
    else 
    return EMPTY;
}

int within_bounds(float value, float min, float max) {
    return value >= min && value <= max;
}
voxel get_cuboid_voxel(int x, int y, int z){
    rave_vec3 cuboid_center = (rave_vec3){-14, 0, -10};
    rave_vec3 cuboid_half_size = (rave_vec3){10, 10, 20};  // half the dimensions of the cuboid
    rave_vec3 vox_center = add3((rave_vec3){x, y, z}, (rave_vec3){0.5f, 0.5f, 0.5f});

    // Check if the voxel's center is within the cuboid's bounds
    if (within_bounds(vox_center.x, cuboid_center.x - cuboid_half_size.x, cuboid_center.x + cuboid_half_size.x) &&
        within_bounds(vox_center.y, cuboid_center.y - cuboid_half_size.y, cuboid_center.y + cuboid_half_size.y) &&
        within_bounds(vox_center.z, cuboid_center.z - cuboid_half_size.z, cuboid_center.z + cuboid_half_size.z)) {
        return BOX;
    } else 
    return EMPTY;
}
voxel get_sphere_voxel(int x, int y, int z){
    rave_vec3 sphere_center = (rave_vec3) {+10,0,-15};
    rave_vec3 vox_center = add3((rave_vec3){x,y,z}, (rave_vec3){0.5,0.5,0.5});
    float d = distance3(vox_center, sphere_center);
    // return 0;
    if(d < 10.0){
        return SPHERE;
    } else return EMPTY;
}
voxel get_light_voxel(int x, int y, int z){
    if(
        (z==29) 
        && 
        (
            ((x>-15) && (x<+15)) 
            && 
            ((y>-15) && (y<+15))
        )) return LIGHT;
    else return EMPTY;
}
voxel rave_get_voxel(int x, int y, int z){
    // rave_vec3 vox_center = (rave_vec3){x,y,z};
    voxel box_voxel = get_bounding_box_voxel(x,y,z);
    voxel sphere_voxel = get_sphere_voxel(x,y,z);
    voxel cuboid_voxel = get_cuboid_voxel(x,y,z);
    voxel light_voxel = get_light_voxel(x,y,z);
    // printf("%d ", box_voxel + sphere_voxel + light_voxel);
    return box_voxel + sphere_voxel + light_voxel + cuboid_voxel;
}

static int width = 512, height = 512;
uint32_t max_reflections = 7;
uint32_t max_steps = 256;
const int sample_count = 100;

static rave_vec3 camera_ray_dir;
static rave_vec3 camera_ray_dir_plane;
static rave_vec3 horizline;
static rave_vec3 vertiline;

rave_vec3 rave_get_ray_pos(int x, int y, int z){
    rave_vec3 pos = {0.0*20.0, 1.0*20.0, 0.2*20.0};

    float uv_x = (x*2 - width ) / (float)width ;
    float uv_y = (y*2 - height) / (float)height;

    float scale_x = uv_x * 20.0;
    float scale_y = uv_y * 20.0;
    
    pos = add3(pos, mul3(horizline, scale_x));
    pos = add3(pos, mul3(vertiline, scale_y));

    return pos;
}
rave_vec3 rave_get_ray_dir(int x, int y, int z){
    float uv_x = (x*2 - width ) / (float)width ;
    float uv_y = (y*2 - height) / (float)height;


    return normalize3(
        add3(
            add3(
                camera_ray_dir,
                mul3(vertiline, uv_y*0.5)
            ),
            mul3(horizline, uv_x*0.5)
        ));
}
Material rave_get_material(voxel voxel){
    Material mat = {0};
    mat.roughness = 0.9;
    mat.emmitance = 0.1;
    switch (voxel) {
        case RED_WALL: {
            mat.color = (rave_vec3){0.9,.1,.1};
            break;
        }
        case GREEN_WALL:{
            mat.color = (rave_vec3){0.05,0.9,0.05};
            break;
        }
        case WHITE_WALL:{
            mat.color = (rave_vec3){0.8,0.8,0.8};
            // mat.roughness = 0.15;
            break;
        }
        case GLOSSY_WALL:{
            mat.color = (rave_vec3){0.9,0.9,0.9};
            break;
        }
        case LIGHT:{
            mat.color = (rave_vec3){0.9,0.9,0.9};
            mat.emmitance = 0.9;
            break;
        }
        case SPHERE:{
            mat.color = (rave_vec3){0.6,0.6,0.6};
            break;
        }
        case BOX:{
            mat.color = (rave_vec3){0.5,0.5,0.5};
            break;
        }
        
    }
    // mat.emmitance = 0.01;
    // if (voxel == 2) {
    //     mat.color = (rave_vec3){0.1,.2,.9};
    //     mat.emmitance = 0.9;
    // }
    // if (voxel == 3) {
    //     mat.color = (rave_vec3){0.1,.9,.1};
    //     mat.emmitance = 0.3;
    // }
    return mat;
}

typedef struct {
    float r,g,b;
} stored_light;
stored_light* image;

static uint8_t *data;
static tga_info *info;
void rave_store_light(int x, int y, int z, rave_vec3 light){

    image[x + width*y].r += (light.x / ((float)sample_count));
    image[x + width*y].g += (light.y / ((float)sample_count));
    image[x + width*y].b += (light.z / ((float)sample_count));
}

int main(int argc, char *argv[]) {
    const char* out_name = "image.tga";
    image = malloc(sizeof(stored_light) * width*height);
    assert(image != NULL);

    //lol never use exactly 0.0
    camera_ray_dir = normalize3((rave_vec3){-0.001, -1.0, -0.2});
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

    for(int w=0; w<width; w++){
    for(int h=0; h<height; h++){
        uint8_t *pixel;
        rave_vec3 light;
        light.x = image[w + width*h].r;
        light.y = image[w + width*h].g;
        light.z = image[w + width*h].b;

        rave_vec3 final_light = clamp3(light, 0, 1);
        pixel = tga_get_pixel(data, info, w, h);
        pixel[0] += (uint8_t)(final_light.x * 255.0);
        pixel[1] += (uint8_t)(final_light.y * 255.0);
        pixel[2] += (uint8_t)(final_light.z * 255.0);
    }}

    // Saves the image as a TGA file.
    error_code = tga_save_from_info(data, info, out_name);
    assert(error_code == TGA_NO_ERROR);

    tga_free_data(data);
    tga_free_info(info);
    free(image);

    return 0;
}