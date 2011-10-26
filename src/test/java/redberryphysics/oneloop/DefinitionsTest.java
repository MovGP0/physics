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
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package redberryphysics.oneloop;

import redberry.core.tensor.Tensor;
import java.util.Deque;
import redberry.core.transformation.collect.CollectPatternParser;
import org.junit.Test;
import redberry.core.tensor.MultiTensor;
import redberry.core.tensor.test.TTest;
import redberry.core.transformation.Transformations;
import redberry.core.transformation.collect.CollectManager;
import redberry.core.transformation.collect.SplitPattern;
import static redberryphysics.core.oneloop.Definitions.*;
import static core.TAssert.*;
import static org.junit.Assert.*;

/**
 *
 * @author stas
 */
public class DefinitionsTest {
    @Test
    public void substitutions1() {
        RR = RR_SUBSTITUTION.transform(RR);
        RR = DELTA_1_SUBSTITUTION.transform(RR);
        RR = DELTA_2_SUBSTITUTION.transform(RR);
        RR = DELTA_3_SUBSTITUTION.transform(RR);
        RR = DELTA_4_SUBSTITUTION.transform(RR);
        RR = Transformations.multiplyNumbers(RR);
        RR = Transformations.sumNumbers(RR);
        Tensor copy = RR.clone();
        assertTrue(TTest.testIsScalar(RR));
        assertIndexes(RR);
        System.out.println(((MultiTensor) RR).size());
        RR = Transformations.expandBrackets(RR);
        assertTrue(TTest.testIsScalar(RR));
        assertIndexes(RR);
        System.out.println(((MultiTensor) RR).size());
        Deque<SplitPattern> deque = CollectPatternParser.parse("(HATK_{\\mu}"
                + "|HATK_{\\mu\\nu}"
                + "|HATK_{\\mu\\nu\\alpha}"
                + "|HATK_{\\mu\\nu\\alpha\\beta}"
                + "|n_{\\mu}"
                + "|(R_{\\mu\\nu}|R_{\\mu\\nu\\alpha\\beta}))");
        CollectManager collectManager = new CollectManager(deque);

        RR = collectManager.collect(RR);
        assertTrue(TTest.testIsScalar(RR));
        assertIndexes(RR);
        System.out.println(((MultiTensor) RR).size());
        RR = Transformations.multiplyNumbers(RR);
        RR = Transformations.sumNumbers(RR);
    }
}
