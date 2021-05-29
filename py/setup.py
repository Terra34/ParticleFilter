from distutils.core import setup, Extension


def main():
    setup(name="generators",
          version="1.0.0",
          description="Generator module with gaussian random generators - Contains: Inverse, Box-Muller, Marsaglia Polar, Ziggurat and Linear feedback shift register implementations.",
          author="<Terra34>",
          author_email="I.Mareticcccc@gmail.com",
          ext_modules=[Extension("generators", ["generator_ext.c"])])


if __name__ == "__main__":
    main()