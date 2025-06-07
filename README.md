# CSE305Project

Final project for CSE 305.

## 🛠️ Getting Started

Clone the repository and initialize its submodules:

```bash
git clone https://github.com/ThisIsErik/CSE305Project.git
cd CSE305Project
git submodule update --init --recursive
```

## 🧱 Build Options
There are two build targets available:

make – Builds the default version.

make cpu – Builds the CPU-only version.

To compile, run:
```bash
make        # Default build
# or
make cpu    # CPU-only build
```
## 🚀 Running the Project
After building, run the executable with a specific mode:

```bash
./run <mode>         # For default build
./run_cpu <mode>     # For CPU-only build
```

Available Modes:

correctness – Runs basic correctness checks for all algorithms.

thread – Benchmarks performance as a function of thread count.

size – Measures runtime and speedup as input size increases.

similarity – Compares performance for different sequence similarities.

memory – Measures memory usage (VmPeak and VmHWM) for selected algorithms.
