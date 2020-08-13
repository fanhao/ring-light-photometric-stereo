# ring-light-photometric-stereo
Classic photometric stereo uses ring light sources.

We present a new method for correcting the bias in traditional photometric stereo when using a ring-light imaging device. 

(1) We analyze the cause of the bias that exists in reconstruction and build surface gradient and height error models, which are prone to be an approximate linear gradient deviation and a related quadratic height deviation. 

(2) We then propose the ring-light photometric stereo compensation method, that uses sparse known 3D points to fit and correct the height deviation for a more accurate 3D heights of the target.

Function information:

(1) 'PS_test6.m' is our main function.

(2) 'PhotometricStereo.m' computes the normal by the traditional photometric stereo method.

(3) 'poisson_solver_function_neumann.m' integrates the gradients to heights.

Reference:

(1) Light Attenuation Model: https://docs.blender.org/manual/en/2.79/render/blender_render/lighting/lights/attenuation.html

(2) Paper&Code: What is the range of surface reconstructions from a gradient field?
