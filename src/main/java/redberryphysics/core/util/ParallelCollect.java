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
package redberryphysics.core.util;

import redberry.core.tensor.Sum;
import redberry.core.tensor.Tensor;
import redberry.core.transformation.Transformation;
import redberry.core.transformation.collect.CollectFactory;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class ParallelCollect implements Transformation {
    public static final ParallelCollect INSTANCE = new ParallelCollect();

    private ParallelCollect() {
    }

    @Override
    public Tensor transform(Tensor tensor) {
        if (!(tensor instanceof Sum))
            return tensor;
        Sum sum = (Sum) tensor;
        int size = sum.size();
        if (size < 10000)
            return CollectFactory.createCollectEqualTerms().transform(sum);
        int subSize = size / 2;
        Sum subSum1 = new Sum();
        Sum subSum2 = new Sum();
        int count = 0;
        for (Tensor t : sum) {
            if (count < subSize)
                subSum1.add(t);
            else
                subSum2.add(t);
            count++;
        }
        Thread[] threads = new Thread[2];

        Sum[] parts = new Sum[]{subSum1, subSum2};
        for (int i = 0; i < 2; i++) {
            Thread t = new ParallelCollectWorker(parts[i], i);
            t.start();
            threads[i] = t;
        }

        for (int i = 0; i < 2; i++)
            try {
                threads[i].join();
            } catch (InterruptedException ex) {
            }

        Sum fin = new Sum();
        fin.add(subSum1);
        fin.add(subSum2);
        return CollectFactory.createCollectEqualTerms().transform(fin);
    }

    private static class ParallelCollectWorker extends Thread {
        private int index;
        private Tensor toCollect;

        public ParallelCollectWorker(Tensor toCollect, int index) {
            this.index = index;
            this.toCollect = toCollect;
        }

        @Override
        public void run() {
            toCollect = CollectFactory.createCollectEqualTerms().transform(toCollect);
        }
    }
}