halmd_add_library(halmd_mdsim_geometries
  cuboid.cpp
  cylinder.cpp
  sphere.cpp
)

halmd_add_modules(
  libhalmd_mdsim_geometries_cuboid
  libhalmd_mdsim_geometries_cylinder
  libhalmd_mdsim_geometries_sphere
)
