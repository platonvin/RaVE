//julia set

#define RAVE_CUSTOM_VOXEL_TYPE int
// #define RAVE_NO_SIMD
#include "fractal_fun.c"

#include "../src/raytracer.h"
#include "../dependencies/tgafunc.h"

#include <stdio.h>

static int width = 256, height = 256;
uint32_t max_reflections = 7;
uint32_t max_steps = 512;
const int sample_count = 20;
rave_voxel rave_empty_voxel = 0;

static rave_vec3 camera_ray_dir;
static rave_vec3 camera_ray_dir_plane;
static rave_vec3 horizline;
static rave_vec3 vertiline;

void delete_if_exists(const char *filename) {
    assert(filename != NULL);
    FILE *file = fopen(filename, "r");
    assert((file && !fclose(file) && !remove(filename)) || !file);
}
//you probably want to have bounding box to prevent long traversals
//you may also enable (uncomment) return on low accumulated reflection and make bounding box black
//also, prefer precomputing voxel grids so it is just 3d array element access
const int BOUNDS = 60;
rave_voxel get_bounding_box_voxel(int x, int y, int z){
    if((x==-BOUNDS) || (y==-BOUNDS) || (z==-BOUNDS) ||
       (x==+BOUNDS) || (y==+BOUNDS) || (z==+BOUNDS)
    ) return 5;
    
    return 0;
}
rave_voxel get_light_voxel(int x, int y, int z){
    if(z== +(BOUNDS-1)) return 16; //just random special value for light
    if(z== -(BOUNDS-1)) return 16; //just random special value for light
    // if( (abs(x) < 2) && 
    //     (abs(y) < 2) && 
    //     (abs(z) < 2)) return 16; //just random special value for light

    return 0;
}
rave_voxel rave_get_voxel(int x, int y, int z){
    //from https://platonvin.github.io/
    // C_x = -0.65 + -0.65*0.15 * Math.sin((currentTime + 23.0)/2.12);
    // C_y = -0.3  + -0.3 *0.15 * Math.sin((currentTime + 23.0)/3.523);
    // C_z = +0.6  + +0.6 *0.15 * Math.sin((currentTime + 23.0)/5.634);
    // C_w = -0.2  + -0.2 *0.15 * Math.sin((currentTime + 23.0)/7.6345);
    
    rave_vec4 Constant = {
        -0.65 + -0.65*0.15 * (sinf(15.0 + 23.0)/2.12  ),
        -0.3  + -0.3 *0.15 * (sinf(15.0 + 23.0)/3.523 ),
        +0.6  + +0.6 *0.15 * (sinf(15.0 + 23.0)/5.634 ),
        -0.2  + -0.2 *0.15 * (sinf(15.0 + 23.0)/7.6345),
        };
    //scaling required, otherwise whole fractal fits into single voxel
    rave_vec3 pos = {
        (((float)x) + 0.5f) / 30.0,
        (((float)y) + 0.5f) / 30.0,
        (((float)z) + 0.5f) / 30.0};

    //uncomment to remove upper half of the fractal
    // if (z < 0) {
        int converged_i = 0;
        float julia = check_julia(pos, Constant, &converged_i);
        if (julia < ESCAPE_THRESHOLD) {
            // printf("js: %f %d\n", julia, converged_i);
            // float division = (0.5+sin(julia * 5624536.0)*0.5)*10.0;
            // return (int)division + 1;
            return (converged_i+3);
        }
    // }

    return get_bounding_box_voxel(x,y,z) + get_light_voxel(x,y,z);
}

rave_vec3 rave_get_ray_pos(int x, int y, int z){
    rave_vec3 pos = rave_mul3(camera_ray_dir, -30);

    float uv_x = (x*2 - width ) / (float)width ;
    float uv_y = (y*2 - height) / (float)height;

    float scale_x = uv_x * 35.0;
    float scale_y = uv_y * 35.0;
    
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
                rave_mul3(vertiline, uv_y*0.1)
            ),
            rave_mul3(horizline, uv_x*0.1)
        ));
}

#define HEX2RBG(_color) \
(rave_vec3){\
    ((((int) _color            ) % 256) / 256.0),\
    ((((int) _color / 256      ) % 256) / 256.0),\
    ((((int) _color / 256 / 256) % 256) / 256.0),\
}

#define FRACTAL_COLORING(X) \
HEX2RBG(0x##X)

#define PALETTE_SIZE 10
rave_vec3 palette_1 [PALETTE_SIZE] = {
    FRACTAL_COLORING(54cea7),
    FRACTAL_COLORING(2ba4a6),
    FRACTAL_COLORING(0c6987),
    FRACTAL_COLORING(054b84),
    FRACTAL_COLORING(0d2147),
    FRACTAL_COLORING(ffb0bf),
    FRACTAL_COLORING(ff82bd),
    FRACTAL_COLORING(d74ac7),
    FRACTAL_COLORING(a825ba),
    FRACTAL_COLORING(682b9c),
};
rave_vec3 palette_2 [PALETTE_SIZE] = {
    FRACTAL_COLORING(fffbef),
    FRACTAL_COLORING(ffcdcd),
    FRACTAL_COLORING(f6857d),
    FRACTAL_COLORING(bd4b6c),
    FRACTAL_COLORING(852646),
    FRACTAL_COLORING(391f46),
    FRACTAL_COLORING(1bdc8d),
    FRACTAL_COLORING(028fa0),
    FRACTAL_COLORING(32d5ff),
    FRACTAL_COLORING(c4fffb),
};
rave_vec3 get_palette_color(int index){
    // return palette_1[index-1];
    return palette_2[index-1];
}

Material rave_get_material(rave_voxel voxel){
    Material mat = {0};
    mat.roughness = 0.95; // set to 0 for mirror

    // mat.color = (rave_vec3){0.65,.725,.8};
    if(voxel == 16){
        mat.emmitance = 4.0;
        mat.color = get_palette_color(6);
    } else {
        mat.color = get_palette_color(voxel);
        // if(voxel == 6){
        //     mat.emmitance = 1.0;
        // }
    }



    return mat;
}

typedef struct {
    float r,g,b;
} stored_light;
stored_light* image;

static uint8_t *data;
static tga_info *info;
#define VERIFY_COLOR(clr) \
assert(clr != NAN);\
assert(clr != +INFINITY);\
assert(clr != -INFINITY);
void rave_store_light(int x, int y, int z, rave_vec3 light){
    // light = rave_clamp3(light, 0, 1);

    // VERIFY_COLOR(image[x + width*y].r);
    // VERIFY_COLOR(image[x + width*y].g);
    // VERIFY_COLOR(image[x + width*y].b);
    image[x + width*y].r += (light.x / ((float)sample_count));
    image[x + width*y].g += (light.y / ((float)sample_count));
    image[x + width*y].b += (light.z / ((float)sample_count));
}

static uint8_t *data;
static tga_info *info;

int main(int argc, char *argv[]) {
    const char* out_name = "fractal_image.tga";
    image = malloc(sizeof(stored_light) * width*height);

    camera_ray_dir = rave_normalize3((rave_vec3){-0.5, -0.8, -0.2});
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

    return 0;
}