{
  description = "A flake for the trustmap repo";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            pkgs.cargo
            # pkgs.rustc
            pkgs.linuxKernel.packages.linux_zen.perf
            pkgs.cargo-cross
            pkgs.cargo-release
            pkgs.htslib

            # For building
            pkgs.openssl
            pkgs.clang
            pkgs.pkg-config
            pkgs.stdenv.cc.cc.lib

            # For python package
            pkgs.python3Full
            pkgs.maturin
          ];

          shellHook = ''
            export LIBCLANG_PATH="${pkgs.llvmPackages.libclang.lib}/lib"
            export LD_LIBRARY_PATH="${pkgs.stdenv.cc.cc.lib}/lib/:/run/opengl-driver/lib/"
          '';
        };

        packages.default = pkgs.buildRustPackage {
          src = ./.;
          cargoLock = {
            lockFile = ./Cargo.lock;
          };
          buildInputs = [
            pkgs.openssl
            pkgs.pkg-config
          ];
          nativeBuildInputs = [
            pkgs.openssl_3_3
            pkgs.pkg-config
          ];
        };
      }
    );
}
