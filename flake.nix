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
        rEnv = pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            tidyverse
            languageserver
            janitor
            data_table
            Rtsne
            Biostrings
          ];
        };
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            rEnv
            pkgs.cargo
            pkgs.rustc
            pkgs.linuxKernel.packages.linux_zen.perf
            pkgs.cargo-cross
            pkgs.cargo-release
            pkgs.htslib
            pkgs.samtools

            # For building
            pkgs.openssl
            pkgs.clang
            pkgs.pkg-config
            pkgs.stdenv.cc.cc.lib

            # For python package
            pkgs.python3
            pkgs.maturin

            # Python dependencies
            (pkgs.python3.withPackages (
              ps: with ps; [
                pysam
                numpy
              ]
            ))

            # System libraries
            pkgs.zlib
            pkgs.bzip2
            pkgs.xz
            pkgs.libdeflate
          ];

          shellHook = ''
            export LIBCLANG_PATH="${pkgs.llvmPackages.libclang.lib}/lib"
            export LD_LIBRARY_PATH="${pkgs.stdenv.cc.cc.lib}/lib/:${pkgs.zlib}/lib:${pkgs.bzip2.out}/lib:${pkgs.xz.out}/lib:${pkgs.libdeflate}/lib:/run/opengl-driver/lib/"

            mkdir -p $HOME/.R/library
            export R_LIBS_USER=$HOME/.R/library
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
