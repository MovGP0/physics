/*
 *
 * Redberry: symbolic tensor computations library.
 * Copyright (C) 2010-2011  Stanislav Poslavsky <stvlpos@mail.ru>
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
package org.redberry.physics.qgr2;

import org.junit.Test;
import redberry.core.context.CC;
import redberry.core.tensor.Tensor;
import static org.redberry.physics.TAssert.*;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class ToFourierTest {
    public ToFourierTest() {
    }

    @Test
    public void test1() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[f[x],x]");
        assertParity(fourier.transform(target), "f[0]");
    }

    @Test
    public void test2() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[a*f[x]*b,x]");
        assertParity(fourier.transform(target), "a*b*f[0]");
    }

    @Test
    public void test3() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[a*f[x]*g[x]*b,x]");
        assertParity(fourier.transform(target), "Integral[a*f[y0]*g[-y0]*b,y0]");
    }

    @Test
    public void test4() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[a*f[x]*g[x]*f[x]*f[x]*X[y1]*b,x]");
        assertParity(fourier.transform(target), "Integral[a*X[y1]*f[y0]*g[y2]*f[y3]*f[-y0-y2-y3]*b,y0,y2,y3]");
    }

    @Test(expected = UnsupportedOperationException.class)
    public void test5() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[a*f[x,y]*g[x]*f[x]*f[x]*X[y1]*b,x]");
        System.out.println(fourier.transform(target));
    }

    /**
     * Derivatives
     */
    @Test
    public void testD1() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[D[f[x],x],x]");
        assertParity(fourier.transform(target), "0");
    }
    
    @Test
    public void testD2() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[f[x]*D[f[x],x],x]");
        System.out.println(fourier.transform(target));
        assertParity(fourier.transform(target), "Integral[f[y0]*I*(-1)*y0*f[(-1)*y0],y0]");
    }
}
