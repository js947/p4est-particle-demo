#include <cstdio>
#include <random>

int main(int argc, char **argv) {
  size_t n = 1000;
  if (argc > 1)
    sscanf(argv[1], "%zd", &n);

  std::random_device r;
  std::default_random_engine gen(r());

  std::normal_distribution<double> x{0.5, 0.1}, y{0.5, 0.1};

  for (size_t i=0; i<n; ++i)
    printf("%f\t%f\n", x(gen), y(gen));
}
