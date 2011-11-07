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
package redberryphysics.core.oneloop.gravity;

import org.junit.Test;
import redberry.core.context.CC;
import redberry.core.tensor.Expression;
import redberry.core.tensor.Tensor;
import redberry.core.tensor.TensorIterator;
import redberry.core.tensor.test.TTest;
import redberry.core.transformation.CalculateNumbers;
import redberry.core.transformation.ExpandBrackets;
import redberry.core.transformation.Transformer;
import redberry.core.transformation.collect.CollectFactory;
import redberry.core.transformation.collect.CollectPowers;
import redberry.core.transformation.contractions.IndexesContractionsTransformation;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class FeynmanRulesTest {
    public FeynmanRulesTest() {
    }

    @Test
    public void test1() {
        FeynmanRules feynmanRules = new FeynmanRules();
        Expression loop1 =
                new Expression("SC[k^a] = T^{mn}[p^a,k^a-p^a]*T_{mn}[p^a-k^a,p^a]");
        loop1.eval(feynmanRules.scalar2_gluon_vertex.asSubstitution());
        loop1.eval(
                new Transformer(ExpandBrackets.EXPAND_ALL),
                IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                feynmanRules.KRONECKER_DIMENSION.asSubstitution(),
                CalculateNumbers.INSTANCE,
                CollectFactory.createCollectAllEqualTerms(),
                CollectFactory.createCollectAllScalars(),
//                new Expression("p_a=x*k_a"),
                CalculateNumbers.INSTANCE,
                new Transformer(CollectPowers.INSTANCE));
        Tensor rhs = loop1.right();
        TensorIterator iterator = rhs.iterator();
        Tensor current;
        int count;
        while (iterator.hasNext()) {
            count = 0;
            current = iterator.next();
            for (Tensor t : current)
                if (TTest.testEqualstensorStructure(t, CC.parse("p^a")))
                    count++;
            if (count % 2 != 0)
                iterator.remove();

        }
        System.out.println(loop1);
    }

    @Test
    public void test2() {
        FeynmanRules feynmanRules = new FeynmanRules();
        Expression e = new Expression("A_mnab = G_mnab[p^q]");
        e.eval(feynmanRules.propagator.asSubstitution(), IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC);
        System.out.println(e);
    }
}
