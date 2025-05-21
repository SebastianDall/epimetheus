{
  description = "A flake for the trustmap repo";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs?ref=nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };
      in {
        devShells.default = pkgs.mkShell {
          buildInputs = [
            pkgs.cargo
            pkgs.rustc
            pkgs.linuxKernel.packages.linux_zen.perf
            pkgs.cargo-cross
            pkgs.cargo-release

            # For building
            pkgs.openssl
            pkgs.clang
            pkgs.pkg-config

            # For interacting with mongodb
            pkgs.mongosh

            # Testing http endpoints
            pkgs.hurl
          ];
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
    });
}
