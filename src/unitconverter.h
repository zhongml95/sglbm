#pragma once
#ifndef UNITCONVERTER_H
#define UNITCONVERTER_H

#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include "parameters.h"

bool readParametersFromXML(const std::string& filePath, Parameters& params);

template<typename T>
class UnitConverter {
public:
  // Public data (initialized by init); you can make these private later if you prefer
  int    N{0};
  int    nx{0};
  T      L{1};
  T      dx{1};
  T      dt{1};
  T      tau{0};

  T      physVelocity{0};
  T      physViscosity{0};
  T      physDensity{1};
  T      Re{0};
  T      u0{0};
  T      omega0{0};

  T      conversionVelocity{1};
  T      conversionViscosity{1};
  T      conversionDensity{1};
  T      conversionMass{1};
  T      conversionForce{1};

  UnitConverter() = default;

  // ----------------- Main entry -----------------
  bool loadFromXML(const std::string& filePath) {
    Parameters params;
    if (!readParametersFromXML(filePath, params)) {
      std::cerr << "[UnitConverter] Failed to load parameters from " << filePath << std::endl;
      return false;
    }
    UnitConverterFromResolutionAndRelaxationTime(params);
    return true;
  }


  // Keep legacy API name
  void UnitConverterFromResolutionAndRelaxationTime(const Parameters& params) {
    N             = static_cast<int>(params.resolution);
    nx            = N;
    // tau           = static_cast<T>(params.tau);
    L             = static_cast<T>(params.L);
    physVelocity  = static_cast<T>(params.physVelocity);
    physDensity   = T(1);

    if (params.Re == 0.0) {
      physViscosity = static_cast<T>(params.physViscosity);
      if (physViscosity == T(0)) {
        throw std::runtime_error("[UnitConverter] Re == 0 and physViscosity == 0: underdetermined.");
      }
      Re = physVelocity * L / physViscosity;
    } else if (params.physViscosity == 0.0) {
      Re = static_cast<T>(params.Re);
      if (Re == T(0)) {
        throw std::runtime_error("[UnitConverter] physViscosity == 0 and Re == 0: underdetermined.");
      }
      physViscosity = physVelocity * L / Re;
    } else {
      physViscosity = static_cast<T>(params.physViscosity);
      Re = physVelocity * L / physViscosity;
    }

    physViscosity = physVelocity * L / Re;
    tau = 3 * physViscosity + 0.5;

    dx = L / static_cast<T>(N);
    dt = (tau - T(0.5)) / T(3.0) * (dx * dx) / physViscosity;

    conversionViscosity = dx * dx / dt;
    conversionVelocity  = dx / dt;
    conversionDensity   = T(1);
    conversionMass      = conversionDensity * dx * dx * dx; // keep dx^3 (3D convention)
    conversionForce     = conversionMass * dx / (dt * dt);
    u0                  = physVelocity / conversionVelocity;
    omega0              = T(1) / (T(3) * (physViscosity / conversionViscosity) + T(0.5));
  }

  void initialize(const Parameters& p) { UnitConverterFromResolutionAndRelaxationTime(p); }

  // ---------- Getter API (const, noexcept) ----------
  // geometry / time
  int    getResolution()           const noexcept { return N; }
  int    getNx()                   const noexcept { return nx; }
  T      getL()                    const noexcept { return L; }
  T      getDx()                   const noexcept { return dx; }
  T      getDt()                   const noexcept { return dt; }
  T      getTau()                  const noexcept { return tau; }

  // physical properties
  T      getPhysVelocity()         const noexcept { return physVelocity; }
  T      getPhysViscosity()        const noexcept { return physViscosity; }
  T      getPhysDensity()          const noexcept { return physDensity; }
  T      getRe()                   const noexcept { return Re; }

  // lattice scales
  T      getU0()                   const noexcept { return u0; }
  T      getOmega0()               const noexcept { return omega0; }
  T      getMa()                   const noexcept { return u0 * std::sqrt(T(3)); } // c_s=1/sqrt(3)

  // conversions
  T      getConversionVelocity()   const noexcept { return conversionVelocity; }
  T      getConversionViscosity()  const noexcept { return conversionViscosity; }
  T      getConversionDensity()    const noexcept { return conversionDensity; }
  T      getConversionMass()       const noexcept { return conversionMass; }
  T      getConversionForce()      const noexcept { return conversionForce; }

  // pretty log (optional)
  void print(std::ostream& os = std::cout, const std::string& title = "Unit Converter") const {
    os << "==== " << title << " ====\n"
       << std::setprecision(12)
       << "N (resolution)         = " << N << "\n"
       << "L                      = " << L << "\n"
       << "tau                    = " << tau << "\n"
       << "dx                     = " << dx << "\n"
       << "dt                     = " << dt << "\n"
       << "physVelocity           = " << physVelocity << "\n"
       << "physViscosity (nu)     = " << physViscosity << "\n"
       << "physDensity            = " << physDensity << "\n"
       << "Re                     = " << Re << "\n"
       << "u0 (lattice)           = " << u0 << "\n"
       << "omega0                 = " << omega0 << "\n"
       << "conversionVelocity     = " << conversionVelocity << "\n"
       << "conversionViscosity    = " << conversionViscosity << "\n"
       << "conversionDensity      = " << conversionDensity << "\n"
       << "conversionMass         = " << conversionMass << "\n"
       << "conversionForce        = " << conversionForce << "\n"
       << "Ma (computed)          = " << getMa() << "\n"
       << "========================\n";
  }
};

#endif // UNITCONVERTER_H
