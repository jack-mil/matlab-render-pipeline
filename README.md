# 3D Rendering Pipeline in MATLAB

## Features
- Perspective or Orthographic projection (Right Handed Axis)
- Full object transform (Scale, Rotation, Translation)
- Freely position camera in 3D space
- Configurable camera FOV and aspect ratio
- Any number of point lights, with any color
- Per vertex shading with Specular & Diffuse lighting passes (ambient optional)
- Only 30 lines of code! Fully vectorized, no loops used

## Images

![rotating suzanne gif](images/spin.gif)

![specular & diffuse shading](images/sphere_phong.png)

![multiple lights](images/sphere_studio.png)

![hi-poly model](images/suzanne_studio.png)

![colored backlighting](images/suzanne_backlit.png)

![utah teapot](images/teapot_studio.png)

![bevelled cube](images/cube_backlit.png)

![suzanne orthographic](images/suzanne_ortho_colored_lights.png)

![orange torus](images/torus_orange.png)

![grey torus](images/torus.png)


Animation made by rendering frames and encoding in various formats with ffmpeg and gifksi.
All actual video formats compress better then the gif, even when lossless. I wish GitHub 
allowed auto-playback of webm/vp9 in embedded markdown.
```bash
gifski -r 30 --extra -Q100 -W728 --repeat 4 -o spin.gif frames/*.png

ffmpeg -framerate 30 -i 'frames/frame%03d.png' -pix_fmt yuv420p -c:v libx264 -preset veryslow -crf 18 -movflags +faststart -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" spin.mp4

ffmpeg -framerate 30 -i frames/frame%03d.png -pix_fmt yuv420p -c:v libvpx-vp9 -b:v 0 -crf 31 -pass 1 -f null /dev/null && \
ffmpeg -framerate 30 -i frames/frame%03d.png -pix_fmt yuv420p -c:v libvpx-vp9 -b:v 0 -crf 31 -pass 2 images/spin-2pass.webm

ffmpeg -framerate 30 -i frames/frame%03d.png -pix_fmt yuv420p -c:v libvpx-vp9 -lossless 1 images/spin-lossless.webm
```

Utah Teapot from ![Wikipedia, under CC0 Public Domain](https://wikipedia.org/wiki/File%3AUtah_teapot_%28solid%29.stl)

All other models created by me in Blender.
