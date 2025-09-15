/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Mingliang Zhong, Stephan Simonis
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

// fft.h

#ifndef FFT_H
#define FFT_H

#include "matrixOperation.h"
#include <cassert>

using Scalar  = double;
using Complex = std::complex<Scalar>;
constexpr Scalar PI = 3.14159265358979323846;


namespace olb {

namespace uq {

namespace fft {

template <typename T>
void fft_dispatch(std::vector<std::complex<T>>& a);

template <typename T>
void ifft_dispatch(std::vector<std::complex<T>>& a);


// Recursive radix-2 Cooley–Tukey FFT
template <typename T>
void fft(std::vector<std::complex<T>>& a) {
    const std::size_t N = a.size();
    if (N <= 1) return;

    if ((N & (N - 1)) != 0) {
        throw std::runtime_error("fft() expects power-of-two size. Use fft_dispatch instead.");
    }

    std::vector<std::complex<T>> even(N / 2), odd(N / 2);
    for (std::size_t i = 0; i < N / 2; ++i) {
        even[i] = a[i * 2];
        odd[i] = a[i * 2 + 1];
    }

    fft(even);
    fft(odd);

    for (std::size_t k = 0; k < N / 2; ++k) {
        std::complex<T> t = std::polar<T>(1.0, -2.0 * M_PI * k / N) * odd[k];
        a[k] = even[k] + t;
        a[k + N / 2] = even[k] - t;
    }
}

// General FFT via Bluestein's algorithm
template <typename T>
void fft_bluestein(std::vector<std::complex<T>>& a) {
    std::size_t N = a.size();
    std::size_t M = 1;
    while (M < 2 * N - 1) {
        M <<= 1;
    }

    std::vector<std::complex<T>> expTable(N);
    for (std::size_t k = 0; k < N; ++k) {
        T angle = M_PI * k * k / N;
        expTable[k] = std::polar<T>(1.0, -angle);
    }

    std::vector<std::complex<T>> A(M, 0.0), B(M, 0.0);
    for (std::size_t k = 0; k < N; ++k) {
        A[k] = a[k] * expTable[k];
        B[k] = std::conj(expTable[k]);
    }
    for (std::size_t k = 1; k < N; ++k) {
        B[M - k] = std::conj(expTable[k]);
    }

    fft_dispatch(A);
    fft_dispatch(B);

    for (std::size_t i = 0; i < M; ++i) {
        A[i] *= B[i];
    }

    ifft_dispatch(A);

    for (std::size_t k = 0; k < N; ++k) {
        a[k] = A[k] * expTable[k];
    }
}

// Dispatch to radix-2 or Bluestein
template <typename T>
void fft_dispatch(std::vector<std::complex<T>>& a) {
    std::size_t N = a.size();
    if ((N & (N - 1)) == 0) {
        fft(a);
    } else {
        fft_bluestein(a);
    }
}


// IFFT (calls fft_dispatch + conjugate normalization)
template <typename T>
void ifft_dispatch(std::vector<std::complex<T>>& a) {
    for (auto& x : a) x = std::conj(x);
    fft_dispatch(a);
    T invN = static_cast<T>(1.0 / a.size());
    for (auto& x : a) x = std::conj(x) * invN;
}

// Real-valued inverse Hermitian FFT (ihfft)
template <typename T>
std::vector<T> ihfft_real(const std::vector<T>& input) {
    std::size_t N = input.size();
    if (N == 0) return {};

    std::size_t M = 2 * N - 2;  // Match numpy.fft.ihfft behavior
    std::vector<std::complex<T>> fullSpectrum(M, std::complex<T>(0.0, 0.0));

    // Fill positive frequency part
    for (std::size_t i = 0; i < N; ++i) {
        fullSpectrum[i] = std::complex<T>(input[i], 0.0);
    }

    // Mirror (conjugate symmetric)
    for (std::size_t i = 1; i < N - 1; ++i) {
        fullSpectrum[M - i] = std::complex<T>(input[i], 0.0);
    }

    // IFFT
    ifft_dispatch(fullSpectrum);

    // Return only the first N values
    std::vector<T> result(N);
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = fullSpectrum[i].real();
    }

    return result;
}



} // namespace FFT

} // namespace uq
} // namespace olb

#endif // FFT_H