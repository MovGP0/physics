/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2012:
 *   Stanislav Poslavsky   <stvlpos@mail.ru>
 *   Bolotin Dmitriy       <bolotin.dmitriy@gmail.com>
 *
 * This file is part of Redberry.
 *
 * Redberry is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Redberry is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Redberry. If not, see <http://www.gnu.org/licenses/>.
 */
package cc.redberry.physics.kv;

import cc.redberry.core.context.CC;
import cc.redberry.core.context.ToStringMode;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.iterators.TensorFirstTreeIterator;
import cc.redberry.core.tensor.testing.TTest;
import cc.redberry.core.utils.Indicator;
import cc.redberry.transformation.CalculateNumbers;
import cc.redberry.transformation.Transformation;
import cc.redberry.transformation.Transformations;
import cc.redberry.transformation.Transformer;
import cc.redberry.transformation.collect.CollecctEqualsInputPort;
import cc.redberry.transformation.collect.CollectFactory;
import cc.redberry.transformation.collect.CollectPowers;
import cc.redberry.transformation.collect.ScalarsSplitCriteria;
import cc.redberry.transformation.concurrent.ExpandAndCollectTransformation;
import cc.redberry.transformation.contractions.IndicesContractionsTransformation;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import cc.redberry.physics.util.SqrSubs;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class MainScenario {
    private static final Transformation EAC = new ExpandAndCollectTransformation(
            new CollecctEqualsInputPort(),
            Indicator.SYMBOL_INDICATOR,
            new Transformation[]{
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                OneLoop.KRONECKER_DIMENSION.asSubstitution(),
                new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}")),
                CalculateNumbers.INSTANCE});
    private static final Transformation SEAC = new ExpandAndCollectTransformation(ScalarsSplitCriteria.INSTANCE,
            Indicator.FALSE_INDICATOR, new Transformation[]{CalculateNumbers.INSTANCE});

    public static void main(String[] args) {
        System.out.println("--------Evaluating " + Integer.parseInt(args[0]) + " term---------");
        OneLoop loop = new OneLoop();
        loop.insertIndices();
        loop.substituteL();
        loop.evalHatK();
        System.out.println(
                ((Sum) loop.RR.right()).getElements().size());
        Delta_Prep.go(loop);
        System.out.println("Deltas done");
        long start = System.currentTimeMillis();
//        for (int i = 0; i < ((Sum) loop.RR.right()).size(); ++i) {
        evalRRTerm(Integer.parseInt(args[0]), loop);
//        }
        long stop = System.currentTimeMillis();
        System.out.println(" TOTAL ---- " + (stop - start));
    }

    /*
     * Calculating fifth element - there is a mistake
     *
     */
//    public static void main(String[] args) {
//        OneLoop loop = new OneLoop();
//        loop.insertIndices();
//        loop.substituteL();
//        loop.evalHatK();
//        System.out.println("Evaluating deltas's");
//        Delta_Prep.go(loop);
//        System.out.println("Done");
//
//        Tensor dSummand = ((Sum) loop.RR.right()).getElements().get(5);
//        long startTime = System.currentTimeMillis();
//        dSummand = dSummand.clone();
//        dSummand = loop.L.asSubstitution().transform(dSummand);
//        dSummand = CalculateNumbers.INSTANCE.transform(dSummand);
//        dSummand = loop.RIMAN.asSubstitution().transform(dSummand);
//        dSummand = loop.RICCI.asSubstitution().transform(dSummand);
//        dSummand = Transformations.expandBrackets(dSummand);
//        dSummand = Transformations.contractMetrics(dSummand);
//        dSummand = CollectFactory.createCollectEqualTerms1().transform(dSummand);
//        dSummand = CalculateNumbers.INSTANCE.transform(dSummand);
//        dSummand = SEAC.transform(dSummand);
//        dSummand = CalculateNumbers.INSTANCE.transform(dSummand);
//        
//        for (Expression h : loop.HATKs)
//            h.asSubstitution().transform(dSummand);
//        for (Expression delta : loop.DELTAs)
//            delta.asSubstitution().transform(dSummand);
//
//        dSummand = Transformations.contractMetrics(dSummand);
//        dSummand = OneLoop.KRONECKER_DIMENSION.asSubstitution().transform(dSummand);
//        dSummand = CalculateNumbers.INSTANCE.transform(dSummand);
//
//        System.out.println("Exp+Collect");
//        dSummand = smartEC(dSummand, EAC);
//
//        System.out.println("+++ Collect Scalar");
//        dSummand = SEAC.transform(dSummand);
//        dSummand = CalculateNumbers.INSTANCE.transform(dSummand);
//        long stopTime = System.currentTimeMillis();
//        System.out.println("+++ Done. Elements: " + getElementsCount(dSummand));
//        System.out.println("+++ Term time = " + (stopTime - startTime) + "ms");
//
//        System.out.println("Evaluating final result ");
//        dSummand = CalculateNumbers.INSTANCE.transform(dSummand);
//        dSummand = SEAC.transform(dSummand);
//        System.out.println("Done: ");
//        System.out.println(dSummand);
//        System.out.println(new Transformer(CollectPowers.INSTANCE).transform(dSummand.clone()));
//
//        if (true)
//            return;
//        System.out.println(CC.parse("L*L*(L-1)*HATK^{\\gamma}*DELTA^{\\alpha\\beta}*HATK^{\\mu\\nu}*n_{\\lambda}*((1/20)*R_{\\alpha\\nu}*R^{\\lambda}_{\\gamma\\beta\\mu}+(1/20)*R_{\\alpha\\gamma}*R^{\\lambda}_{\\mu\\beta\\nu}+(1/10)*R_{\\alpha\\beta}*R^{\\lambda}_{\\mu\\gamma\\nu}+(1/20)*R^{\\sigma}_{\\alpha\\nu\\gamma}*R^{\\lambda}_{\\sigma\\beta\\mu}-(1/60)*R^{\\sigma}_{\\mu\\alpha\\nu}*R^{\\lambda}_{\\beta\\sigma\\gamma}+(1/10)*R^{\\sigma}_{\\alpha\\beta\\gamma}*R^{\\lambda}_{\\mu\\sigma\\nu}-(1/12)*R^{\\sigma}_{\\alpha\\beta\\nu}*R^{\\lambda}_{\\mu\\sigma\\gamma})").toString(ToStringMode.UTF8));
//       
//    }
    private static Tensor evalAll() {

        OneLoop loop = new OneLoop();
        loop.insertIndices();
        loop.substituteL();
        loop.evalHatK();
        System.out.println("Evaluating deltas's");
        Delta_Prep.go(loop);
        System.out.println("Done");

        Sum result = new Sum();
        int RRSize = ((Sum) loop.RR.right()).size();
        long[] times = new long[RRSize];
        int[] elementsCount = new int[RRSize];
        int count = 0;

        for (Tensor dSummand : loop.RR.right()) {
            long startTime = System.currentTimeMillis();
            dSummand = dSummand.clone();
            dSummand = loop.L.asSubstitution().transform(dSummand);
            dSummand = CalculateNumbers.INSTANCE.transform(dSummand);
            dSummand = loop.RIMAN.asSubstitution().transform(dSummand);
            dSummand = loop.RICCI.asSubstitution().transform(dSummand);
            dSummand = Transformations.expandBrackets(dSummand);
            dSummand = Transformations.contractMetrics(dSummand);
            dSummand = CollectFactory.createCollectEqualTerms().transform(dSummand);
            dSummand = CalculateNumbers.INSTANCE.transform(dSummand);

            for (Expression h : loop.HATKs)
                h.asSubstitution().transform(dSummand);
            for (Expression delta : loop.DELTAs)
                delta.asSubstitution().transform(dSummand);

            dSummand = Transformations.contractMetrics(dSummand);
            dSummand = OneLoop.KRONECKER_DIMENSION.asSubstitution().transform(dSummand);
            dSummand = CalculateNumbers.INSTANCE.transform(dSummand);

            System.out.println("Exp+Collect - " + count);
            dSummand = smartEC(dSummand, EAC);

            System.out.println("+++ Collect Scalar - " + count);
            dSummand = SEAC.transform(dSummand);
            dSummand = CalculateNumbers.INSTANCE.transform(dSummand);
            long stopTime = System.currentTimeMillis();
            System.out.println("+++ Done. " + count + " elements: " + getElementsCount(dSummand));
            elementsCount[count] = getElementsCount(dSummand);
            System.out.println("+++ " + count + " term time = " + (stopTime - startTime) + "ms");
            times[count] = stopTime - startTime;
            count++;
            result.add(dSummand);
        }
        System.out.println("Evaluating final result ");
        Tensor finalResult = CalculateNumbers.INSTANCE.transform(result);
        finalResult = SEAC.transform(result);
        System.out.println("Done: ");
        System.out.println(finalResult);
        System.out.println(new Transformer(CollectPowers.INSTANCE).transform(finalResult.clone()));
        try (PrintStream str = new PrintStream("RRresult.txt")) {
            str.println(finalResult);
            str.println("");
            str.println(new Transformer(CollectPowers.INSTANCE).transform(finalResult.clone()));
        } catch (IOException ex) {
        }
        System.out.println("Times : " + Arrays.toString(times));
        System.out.println("Elements : " + Arrays.toString(elementsCount));
        return finalResult;
    }

    static Tensor evalRRTerm(int termIndex, OneLoop loop) {
        Tensor dSummand = ((Sum) loop.RR.right()).getElements().get(termIndex);
        long startTime = System.currentTimeMillis();
        dSummand = dSummand.clone();
        dSummand = loop.L.asSubstitution().transform(dSummand);
        dSummand = CalculateNumbers.INSTANCE.transform(dSummand);
        dSummand = loop.RIMAN.asSubstitution().transform(dSummand);
        dSummand = loop.RICCI.asSubstitution().transform(dSummand);
        dSummand = Transformations.expandBrackets(dSummand);
        dSummand = Transformations.contractMetrics(dSummand);
        dSummand = CollectFactory.createCollectEqualTerms().transform(dSummand);
        dSummand = CalculateNumbers.INSTANCE.transform(dSummand);

        for (Expression h : loop.HATKs)
            h.asSubstitution().transform(dSummand);
        for (Expression delta : loop.DELTAs)
            delta.asSubstitution().transform(dSummand);

        dSummand = Transformations.contractMetrics(dSummand);
        dSummand = OneLoop.KRONECKER_DIMENSION.asSubstitution().transform(dSummand);
        dSummand = CalculateNumbers.INSTANCE.transform(dSummand);

        System.out.println("Exp+Collect");
        dSummand = smartEC(dSummand, EAC);

        System.out.println("+++ Collect Scalar");
        dSummand = SEAC.transform(dSummand);
        dSummand = CalculateNumbers.INSTANCE.transform(dSummand);
        long stopTime = System.currentTimeMillis();
        System.out.println("+++ Done. Elements: " + getElementsCount(dSummand));
        System.out.println("+++ Term time = " + (stopTime - startTime) + "ms");

        System.out.println("Evaluating final result ");
        dSummand = CalculateNumbers.INSTANCE.transform(dSummand);
        dSummand = SEAC.transform(dSummand);
        System.out.println("Done: ");
//        System.out.println(dSummand);
        System.out.println("term = " + new Transformer(CollectPowers.INSTANCE).transform(dSummand.clone()));
        return dSummand;
    }

    static Tensor smartEC(Tensor tensor, Transformation ec) {
        if (tensor instanceof Sum) {
            Sum sum = (Sum) tensor;
            List<Tensor> result = new ArrayList<>();
            int i = 0;
            for (Tensor summand : sum)
//                System.out.println(summand.toString(ToStringMode.UTF8));
                result.add(smartEC(summand, ec)); //                System.out.println(result.get(i++));
            return new Sum(result);
        }
        if (!(tensor instanceof Product))
            throw new IllegalArgumentException();
        ArrayList<Tensor> sums = new ArrayList<>();
        TensorIterator iterator = tensor.iterator();
        while (iterator.hasNext()) {
            Tensor t = iterator.next();
            if (t instanceof Sum && !TTest.testIsSymbol(t)) {
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
                            Transformations.expandAndCollectAllScalars(
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
    private static long fSLEEP_INTERVAL = 100;

    static int getElementsCount(Tensor t) {
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
