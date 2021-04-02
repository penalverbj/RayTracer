# RayTracer
Ray Tracer built as a project for CS 3451 Computer Graphics with Processing p5.js

## Functionality
This ray tracer includes:

1) Intersections with spheres and disks
2) Distributive ray tracing - having several scattered rays shot through the same pixel and averaging the color results for anti-aliasing
3) Shadows - supports both soft and hard shadows. Hard shadows come from point lights and by treating area lights as point lights. Soft shadows come from distributed rays from area lights
4) Jittering - Makes the effect of soft shadows more even by randomizing to what part of the area light the light ray is shot
