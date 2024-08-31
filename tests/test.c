#include "../dependencies/tgafunc.h"

#define RAVE_CUSTOM_VOXEL_TYPE int
#include "../src/raytracer.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

static int width = 256, height = 256;
uint32_t max_reflections = 12;
uint32_t max_steps = 256;
const int sample_count = 10;
rave_voxel rave_empty_voxel = 0;

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
    WHITE_CEILING,
    LIGHT,
    SPHERE,
    BOX,
};
//you probably want to have bounding box to prevent long traversals
//you may also enable (uncomment) return on low accumulated reflection and make bounding box black (material with {0,0,0} color)
//also, prefer precomputing voxel grids so it is just 3d array element access
const int BOUNDS = 30;
rave_voxel get_bounding_box_voxel(int x, int y, int z){
    if(x == -BOUNDS) {
        return RED_WALL;
    } else 
    if (x == +BOUNDS) {
        return GREEN_WALL;
    } else 
    if(z == +BOUNDS){
        return WHITE_CEILING;
    } else
    if ((x <= -BOUNDS) || (y <= -BOUNDS) || (z <= -BOUNDS) || 
        (x >= +BOUNDS) || (y >= +BOUNDS) || (z >= +BOUNDS)) {
        return WHITE_WALL;
    }
    else 
    return EMPTY;
}

int within_bounds(float value, float min, float max) {
    return value >= min && value <= max;
}
rave_voxel get_cuboid_voxel(int x, int y, int z){
    rave_vec3 cuboid_center = (rave_vec3){-14, 0, -10};
    rave_vec3 cuboid_half_size = (rave_vec3){10, 10, 20};
    rave_vec3 vox_center = rave_add3((rave_vec3){x, y, z}, (rave_vec3){0.5f, 0.5f, 0.5f});
    
    if (within_bounds(vox_center.x, cuboid_center.x - cuboid_half_size.x, cuboid_center.x + cuboid_half_size.x) &&
        within_bounds(vox_center.y, cuboid_center.y - cuboid_half_size.y, cuboid_center.y + cuboid_half_size.y) &&
        within_bounds(vox_center.z, cuboid_center.z - cuboid_half_size.z, cuboid_center.z + cuboid_half_size.z)) {
        return BOX;
    } else 
    return EMPTY;
}
rave_voxel get_sphere_voxel(int x, int y, int z){
    rave_vec3 sphere_center = (rave_vec3) {+10,0,-15};
    rave_vec3 vox_center = rave_add3((rave_vec3){x,y,z}, (rave_vec3){0.5,0.5,0.5});
    float d = rave_distance3(vox_center, sphere_center);
    if(d < 10.0){
        return SPHERE;
    } else return EMPTY;
}
rave_voxel get_light_voxel(int x, int y, int z){
    if(
        (z == (BOUNDS-1)) 
        && 
        (
            ((x > -(BOUNDS/2)) && (x < +(BOUNDS/2))) 
            && 
            ((y > -(BOUNDS/2)) && (y < +(BOUNDS/2)))
        )) return LIGHT;
    else return EMPTY;
}
rave_voxel rave_get_voxel(int x, int y, int z){
    rave_voxel box_voxel = get_bounding_box_voxel(x,y,z);
    rave_voxel sphere_voxel = get_sphere_voxel(x,y,z);
    rave_voxel cuboid_voxel = get_cuboid_voxel(x,y,z);
    rave_voxel light_voxel = get_light_voxel(x,y,z);
    return box_voxel + sphere_voxel + light_voxel + cuboid_voxel;
}

static rave_vec3 camera_ray_dir;
static rave_vec3 camera_ray_dir_plane;
static rave_vec3 horizline;
static rave_vec3 vertiline;

rave_vec3 rave_get_ray_pos(int x, int y, int z){
    rave_vec3 pos = {0.001*20.0, 1.0*20.0, 0.2*20.0};

    float uv_x = (x*2 - width ) / (float)width ;
    float uv_y = (y*2 - height) / (float)height;

    float scale_x = uv_x * 20.0;
    float scale_y = uv_y * 20.0;
    
    pos = rave_add3(pos, rave_mul3(horizline, scale_x));
    pos = rave_add3(pos, rave_mul3(vertiline, scale_y));

    return pos;
}
rave_vec3 rave_get_ray_dir(int x, int y, int z){
    float uv_x = (x*2 - width ) / (float)width ;
    float uv_y = (y*2 - height) / (float)height;


    return rave_normalize3(
        rave_add3(
            rave_add3(
                camera_ray_dir,
                rave_mul3(vertiline, uv_y*0.5)
            ),
            rave_mul3(horizline, uv_x*0.5)
        ));
}
Material rave_get_material(rave_voxel voxel){
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
            // mat.roughness = 0.00;
            break;
        }
        case WHITE_CEILING:{
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
    const char* out_name = "Cornell_box_image.tga";
    image = malloc(sizeof(stored_light) * width*height);
    assert(image != NULL);

    //lol never use exactly 0.0
    camera_ray_dir = rave_normalize3((rave_vec3){-0.001, -1.0, -0.2});
    camera_ray_dir_plane = rave_normalize3((rave_vec3){camera_ray_dir.x, camera_ray_dir.y, 0});
    horizline = rave_normalize3(rave_cross3(camera_ray_dir_plane, (rave_vec3){0,0,1}));
    vertiline = rave_normalize3(rave_cross3(camera_ray_dir, horizline));

    delete_if_exists(out_name);
    enum tga_error error_code;
    error_code = tga_create(&data, &info, width, height, TGA_PIXEL_RGB24);
    assert(error_code == TGA_NO_ERROR);

    rave_dispatch(width, height, sample_count);    

    for(int w=0; w<width; w++){
    for(int h=0; h<height; h++){
        uint8_t *pixel;
        rave_vec3 light;
        light.x = image[w + width*h].r;
        light.y = image[w + width*h].g;
        light.z = image[w + width*h].b;

        rave_vec3 final_light = rave_clamp3(light, 0, 1);
        pixel = tga_get_pixel(data, info, w, h);
        pixel[0] += (uint8_t)(final_light.x * 255.0);
        pixel[1] += (uint8_t)(final_light.y * 255.0);
        pixel[2] += (uint8_t)(final_light.z * 255.0);
    }}

    error_code = tga_save_from_info(data, info, out_name);
    assert(error_code == TGA_NO_ERROR);
    tga_free_data(data);
    tga_free_info(info);
    free(image);

    return 0;
}