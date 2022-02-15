# DHB

This is the repository of the Dynamic Hashed Blocks (DHB) format.

## Authors

| Name                      | E-Mail                                  | Affiliation |
|---------------------------|-----------------------------------------|-------------|
| Alexander van der Grinten | avdgrinten@hu-berlin.de                 | HU Berlin   |
| Maria Predari             | predarim@informatik.hu-berlin.de        | HU Berlin   |
| Florian Willich           | florian.willich@informatik.hu-berlin.de | HU Berlin   |

# Build Library

Build the library using CMake (e.g., with the generator Ninja).

```bash
mkdir build
cd build
cmake -GNinja -DCMAKE_BUILD_TYPE=Release ..
ninja
```

# Build Tests

You can simply build the tests by setting the CMake option `DHB_TEST` to `On`.

```
cmake -GNinja -DCMAKE_BUILD_TYPE=Release -DDHB_TEST=On ..
```

# Benchmarks and Experimental Evaluation

An evaluation of DHB and other dynamic graph data structures can be found at: https://github.com/hu-macsy/dhb-experiments
