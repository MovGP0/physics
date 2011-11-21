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
package redberryphysics.core.oneloop;

import java.util.List;
import redberry.core.tensor.iterators.TensorLastTreeIterator;
import redberry.core.tensor.iterators.TensorFirstTreeIterator;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import redberry.core.tensor.SimpleTensor;
import redberry.core.transformation.collect.EqualsSplitCriteria;
import redberry.core.transformation.contractions.IndexesContractionsTransformation;
import redberry.core.utils.Indicator;
import redberryphysics.core.util.SqrSubs;
import redberry.core.context.CC;
import redberry.core.context.ToStringMode;
import redberry.core.tensor.Expression;
import redberry.core.tensor.Product;
import redberry.core.tensor.Sum;
import redberry.core.tensor.Tensor;
import redberry.core.tensor.TensorIterator;
import redberry.core.tensor.TensorWrapper;
import redberry.core.tensor.iterators.TensorTreeIterator;
import redberry.core.tensor.test.TTest;
import redberry.core.transformation.CalculateNumbers;
import redberry.core.transformation.IndexesInsertion;
import redberry.core.transformation.RenameConflictingIndexes;
import redberry.core.transformation.Transformation;
import redberry.core.transformation.Transformations;
import redberry.core.transformation.Transformer;
import redberry.core.transformation.collect.CollecctEqualsInputPort;
import redberry.core.transformation.collect.CollectFactory;
import redberry.core.transformation.collect.ScalarsSplitCriteria;
import redberry.core.transformation.concurrent.ExpandAndCollectTransformation;
import redberry.core.transformation.substitutions.NaiveSubstitution;
import static redberryphysics.core.util.IndexesFactoryUtil.*;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class MainScenario {
    public static void main(String[] args) {
        OneLoop loop = new OneLoop();

        //Only Delta & HATK
        Transformation indexesInsertion;
        indexesInsertion = new IndexesInsertion(loop.matricesIndicator, createIndexes(loop.DELTAs, "^{\\mu\\nu}_{\\alpha\\beta}"));
        loop.DELTA_4.eval(
                indexesInsertion,
                loop.L.asSubstitution(),
                CalculateNumbers.INSTANCE,
                //new Transformer(ExpandBrackets.EXPAND_ALL),
                new Transformer(RenameConflictingIndexes.INSTANCE));

        loop.evalHatK();

        for (Expression ex : loop.HATKs)
            ex.eval(CollectScalars.INSTANCE);

        for (Tensor t : loop.HATKs)
            System.out.println(t.toString(ToStringMode.UTF8));

        Tensor firstSummand = ((Sum) loop.RR.right()).getElements().get(0);

        System.out.println(firstSummand.toString(ToStringMode.UTF8));
        indexesInsertion = new IndexesInsertion(loop.matricesIndicator, doubleAndDumpIndexes(createIndexes(loop.TERMs, "^{\\mu\\nu}")));
        firstSummand = indexesInsertion.transform(firstSummand);
        firstSummand = loop.L.asSubstitution().transform(firstSummand);
        firstSummand = CalculateNumbers.INSTANCE.transform(firstSummand);
        firstSummand = loop.RIMAN.asSubstitution().transform(firstSummand);
        firstSummand = loop.RICCI.asSubstitution().transform(firstSummand);
        firstSummand = Transformations.expandBrackets(firstSummand);
        firstSummand = Transformations.contractMetrics(firstSummand);
        firstSummand = CollectFactory.createCollectEqualTerms1().transform(firstSummand);

        System.out.println(firstSummand.toString(ToStringMode.UTF8));

        Tensor[] hatkCombinations = CC.parse("HATK^{\\alpha\\beta}*HATK^{\\mu\\nu}",
                "HATK^{\\beta}*HATK^{\\alpha}*HATK^{\\mu}*HATK^{\\nu}",
                "HATK^{\\mu}*HATK^{\\nu}*HATK^{\\alpha\\beta}",
                "HATK^{\\mu}*HATK^{\\alpha\\beta}*HATK^{\\nu}",
                "HATK^{\\alpha\\beta}*HATK^{\\mu}*HATK^{\\nu}");
        indexesInsertion = new IndexesInsertion(loop.matricesIndicator, createIndexes(hatkCombinations, "^{\\mu\\nu}_{\\alpha\\beta}"));
        for (int i = 0; i < hatkCombinations.length; ++i) {
            hatkCombinations[i] = indexesInsertion.transform(hatkCombinations[i]);
            System.out.println(hatkCombinations[i].toString(ToStringMode.UTF8));
        }

        Transformation ec = new ExpandAndCollectTransformation(
                new CollecctEqualsInputPort(),
                Indicator.SYMBOL_INDICATOR,
                new Transformation[]{
                    IndexesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    loop.KRONECKER_DIMENSION.asSubstitution(),
                    new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}")),
                    CalculateNumbers.INSTANCE});

        Expression[] hatkCombDone = new Expression[hatkCombinations.length];
        try (PrintStream str = new PrintStream("tensors.txt")) {
            for (int i = 0; i < hatkCombinations.length; ++i) {

                Tensor hatkExpanded = hatkCombinations[i].clone();
                for (Expression hK : loop.HATKs)
                    hK.asSubstitution().transform(hatkExpanded);
                hatkExpanded = smartEC(hatkExpanded, ec);
                hatkExpanded = CollectScalars.INSTANCE.transform(hatkExpanded);
                // hatkExpanded = new Transformer(CollectPowers.INSTANCE).transform(hatkExpanded);
                hatkCombDone[i] = new Expression(hatkCombinations[i], hatkExpanded);
                //if (TensorUtils.testParentConsistent(hatkCombDone[i]))
                //    System.out.println("AAAAAAAAAAAA");
                System.out.println("Mem: " + getMemoryUse());
                str.println(hatkCombDone[i].toString(ToStringMode.UTF8));
            }
        } catch (IOException ex) {
        }

        Tensor cur;
        int count;
        for (Tensor t : hatkCombDone)
            System.out.println("Elements in HATKComb: " + getElementsCount(t));

        System.out.println("Substitution");

        for (Expression t : hatkCombDone)
            NaiveSubstitution.create(t).transform(loop.DELTA_4);

        System.out.println("Mem: " + getMemoryUse());
        System.out.println("Elements in Delta: " + getElementsCount(loop.DELTA_4));

        System.out.println("Collect");

        loop.DELTA_4.eval(ec);

        System.out.println("Mem: " + getMemoryUse());
        System.out.println("Elements in Delta: " + getElementsCount(loop.DELTA_4));

        System.out.println("Collect Scalar");

        CollectScalars.INSTANCE.transform(loop.DELTA_4);

        System.out.println("Mem: " + getMemoryUse());
        System.out.println("Elements in Delta: " + getElementsCount(loop.DELTA_4));

        for (Expression h : loop.HATKs)
            h.asSubstitution().transform(loop.DELTA_4);

        System.out.println("Mem: " + getMemoryUse());
        System.out.println("Elements in Delta: " + getElementsCount(loop.DELTA_4));

        Tensor dc = loop.DELTA_4.clone();
        for (TensorTreeIterator it = new TensorFirstTreeIterator(dc); it.hasNext();)
            if (TTest.testIsSymbol(it.next()))
                it.remove();
        System.out.println(dc.toString(ToStringMode.UTF8));

        System.out.println("Sustitution to FIRST");

        for (Expression h : loop.HATKs)
            h.asSubstitution().transform(firstSummand);

        loop.DELTA_4.asSubstitution().transform(firstSummand);

        firstSummand = Transformations.contractMetrics(firstSummand);
        firstSummand = loop.KRONECKER_DIMENSION.asSubstitution().transform(firstSummand);
        firstSummand = CalculateNumbers.INSTANCE.transform(firstSummand);

        System.out.println(((Sum) loop.DELTA_4.right()).size());
        System.out.println(((Sum) loop.HATK_1.right()).size());

        System.out.println("First elements: " + getElementsCount(firstSummand));

        System.out.println("Exp+Collect");

        firstSummand = smartEC(firstSummand, ec);

        System.out.println(firstSummand);
        System.out.println("Collect Scalar");
        Transformation sc = new ExpandAndCollectTransformation(ScalarsSplitCriteria.INSTANCE,
                Indicator.FALSE_INDICATOR, new Transformation[]{CalculateNumbers.INSTANCE});
        firstSummand = sc.transform(firstSummand);

        System.out.println(firstSummand);
        System.out.println("Done. Elements: " + getElementsCount(firstSummand));
    }
    private static final Transformation scalarsCollectPort =
            new ExpandAndCollectTransformation(
            ScalarsSplitCriteria.INSTANCE,
            Indicator.FALSE_INDICATOR,
            new Transformation[]{CalculateNumbers.INSTANCE});

    private static class CollectScalars implements Transformation {
        public static final CollectScalars INSTANCE = new CollectScalars();

        private CollectScalars() {
        }

        @Override
        public Tensor transform(Tensor tensor) {
            tensor = CalculateNumbers.INSTANCE.transform(tensor);
            Tensor parent = tensor.getParent();
            TensorWrapper wrapper = new TensorWrapper(tensor);
            TensorLastTreeIterator iterator = new TensorLastTreeIterator(wrapper);
            Tensor current;
            while (iterator.hasNext()) {
                current = iterator.next();
                if (TTest.testIsSymbol(current)) {
                    current = scalarsCollectPort.transform(current);
                    current = CalculateNumbers.INSTANCE.transform(current);
                    iterator.set(current);
                }
            }
            tensor = wrapper.getInnerTensor();
            tensor.setParent(parent);
            return tensor;
        }
    }

    private static Tensor smartEC(Tensor tensor, Transformation ec) {
        if (tensor instanceof Sum) {
            Sum sum = (Sum) tensor;
            List<Tensor> result = new ArrayList<>();
            for (Tensor summand : sum)
                result.add(smartEC(summand, ec));
            return new Sum(result);
        }
        if (!(tensor instanceof Product))
            throw new IllegalArgumentException();
        ArrayList<Tensor> sums = new ArrayList<>();
        TensorIterator iterator = tensor.iterator();
        while (iterator.hasNext()) {
            Tensor t = iterator.next();
            if (t instanceof Sum) {
                sums.add(t);
                iterator.remove();
            }
        }
        ArrayList<Tensor> sumsNext = new ArrayList<>(), sumsTmp;
        while (sums.size() > 1) {
            for (int i = 0; i < sums.size() / 2; ++i)
                if (i * 2 + 1 < sums.size()) {
                    System.out.print("Iter: " + ((Sum) sums.get(i * 2)).size() + ", " + ((Sum) sums.get(i * 2 + 1)).size() + "; ");
                    sumsNext.add(
                            CollectScalars.INSTANCE.transform(
                            ec.transform(new Product(sums.get(i * 2), sums.get(i * 2 + 1)))));
                }
            if (sums.size() % 2 == 1)
                sumsNext.add(sums.get(sums.size() - 1));
            sumsTmp = sums;
            sums = sumsNext;
            sumsNext = sumsTmp;
            sumsNext.clear();
            System.out.println("Done.");
        }
        if (((Product) tensor).isEmpty())
            return sums.get(0);
        ((Product) tensor).add(sums.get(0));
        return ec.transform(tensor);
    }
    // PRIVATE //
    private static int fSAMPLE_SIZE = 100;
    private static long fSLEEP_INTERVAL = 100;

    private static int getElementsCount(Tensor t) {
        int count = 0;
        TensorFirstTreeIterator iterator = new TensorFirstTreeIterator(t);
        Tensor cur;
        while (iterator.hasNext()) {
            cur = iterator.next();
//                if (TTest.testIsSymbol(cur))
//                    System.out.println(cur);
//                    iterator.skip();
            count++;
        }
        return count;
    }

    private static long getMemoryUse() {
        putOutTheGarbage();
        long totalMemory = Runtime.getRuntime().totalMemory();

        putOutTheGarbage();
        long freeMemory = Runtime.getRuntime().freeMemory();

        return (totalMemory - freeMemory);
    }

    private static void putOutTheGarbage() {
        collectGarbage();
        collectGarbage();
    }

    private static void collectGarbage() {
        try {
            System.gc();
            Thread.currentThread().sleep(fSLEEP_INTERVAL);
            System.runFinalization();
            Thread.currentThread().sleep(fSLEEP_INTERVAL);
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }
    }
}
