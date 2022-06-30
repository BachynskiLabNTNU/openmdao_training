# openmdao_training

This repo contains slides and tutorial scripts to help users get started with NASA's [OpenMDAO](https://github.com/OpenMDAO/openmdao).  This version is modified from the [original](https://github.com/OpenMDAO/openmdao_training) (with new Lab 0 and Lab 1 examples) for the [SFI BLUES](sfiblues.com) Optimization Workshop, held in June 2022.

When given in a workshop, this material should take about 6 hours to complete.  With this repo, you will be able to go through the slides and scripts at your own pace.

This is not an exhaustive course in design optimization or use of OpenMDAO, and is mainly intended to demonstrate a sample of what is possible with these tools.  Lab 0 and Lab 1 focus on model-building and solvers, while Lab 2 introduces optimization and compares various algorithms.

### To get started:

1. If you don't already have Python installed, download and install [Anaconda Python 3.10](https://www.anaconda.com/distribution/).
2. After that, follow [OpenMDAO's installation instructions](https://github.com/OpenMDAO/openmdao).  Test your installation by running `python paraboloid.py` and `python paraboloid_min.py` in your terminal window (or running each script in an IDE) and comparing with the results below.
3. Then you can go through the `openmdao_workshop_blues` presentation, which will guide you through concepts within OpenMDAO and walk you through the tutorial scripts.
4. Contact Peter Rohrer ([peter.j.rohrer@ntnu.no](mailto:peter.j.rohrer@ntnu.no)) if you have further questions.

#### Check your installation of OpenMDAO

Running `python paraboloid.py` should print the follow to the command prompt:

```
f = -15.00 at x = 3.00, y = -4.00
f = -5.00 at x = 5.00, y = -2.00
```

Running `python paraboloid_min.py` should print the follow to the command prompt:

```
Optimization terminated successfully    (Exit mode 0)
            Current function value: -27.333333333333
            Iterations: 5
            Function evaluations: 6
            Gradient evaluations: 5
Optimization Complete
-----------------------------------
Minimum: -27.33 at x = 6.67, y = -7.33
```

`python paraboloid_min.py` should also generate an N2 diagram, saving a file to the directory (for all setups) and opening the diagram in a browser (if not using an IDE).