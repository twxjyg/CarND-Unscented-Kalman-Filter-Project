# How to coding use VS Code

1. copy .vscode folder into this repo's root directory

    ```shell
        cp .vscode ../../ -r
    ```

2. make build folder to compile and generate compile_commands.json for better c++ auto-completion

    ```shell
        mkdir build
        cd build
        cmake ..
        make
    ```
3. start editting in your VS Code