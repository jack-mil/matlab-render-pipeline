# 3D Rendering Pipeline in MATLAB

### Features
- Perspective or Orthographic projection (Right Handed Axis)
- Full object transform (Scale, Rotation, Translation)
- Freely position camera in 3D space
- Configurable camera FOV and aspect ratio
- Any number of point lights, with any color
- Per vertex shading with Specular & Diffuse lighting passes (ambient optional)
- Only 30 lines of code! Fully vectorized, no loops used

## Images

![specular & diffuse shading](images/sphere_phong.png)
![multiple lights](images/sphere_studio.png)
![hi-poly model](images/suzanne_studio.png)
![colored backlighting](images/suzanne_backlit.png)
![utah teapot](images/teapot_studio.png)
![alt text](images/cube_backlit.png)
![alt text](images/suzanne_ortho_colored_lights.png)
![alt text](images/torus_orange.png)
![alt text](images/torus.png)

![alt text](images/spin.gif)

Animation made by rendering frames and encoding with ffmpeg and gifksi
```bash
gifski -r 30 --extra -Q100 -W728 --repeat 4 -o spin.gif frames/*.png

ffmpeg -framerate 30 -i 'frames/frame%03d.png' -pix_fmt yuv420p -c:v libx264 -preset veryslow -crf 18 -movflags +faststart -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" spin.mp4
```