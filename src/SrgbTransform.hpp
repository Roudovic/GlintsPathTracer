//
//  SrgbTransform.hpp
//  PathTracer
//
//  Created by Ludovic Theobald on 18/12/2019.
//  Copyright © 2019 Ludovic Theobald. All rights reserved.
//

#ifndef SrgbTransform_hpp
#define SrgbTransform_hpp

/*
 * sRGB transform (C++)
 *
 * Copyright (c) 2017 Project Nayuki. (MIT License)
 * https://www.nayuki.io/page/srgb-transform-library
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */

#pragma once


namespace SrgbTransform {

float srgbToLinear(float x);
double srgbToLinear(double x);
extern const float SRGB_8BIT_TO_LINEAR_FLOAT[1 << 8];
extern const double SRGB_8BIT_TO_LINEAR_DOUBLE[1 << 8];

float linearToSrgb(float x);
double linearToSrgb(double x);
int linearToSrgb8bit(double x);

}
#endif /* SrgbTransform_hpp */
