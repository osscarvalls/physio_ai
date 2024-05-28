git clone --recursive -j8 https://github.com/abetlen/llama-cpp-python.git
set FORCE_CMAKE=1
set CMAKE_ARGS=-DLLAMA_CUBLAS=ON
cd llama-cpp-python
python -m pip install -e .