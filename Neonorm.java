/**
 * Used to determine NeONORM bias coefficient for two vectors.
 * 
 * @author Sebastian Noth
 */
public class Neonorm {

	public float lnQ[];
	public float W[];
	private String name;

	public void setName(String n) {
		name = n;
	}

	public Neonorm() {
	}

	public static String getVersion() {
		return "NeONORM 1.0 [050703]";
	}

	private int SGN(double x) {
		return x == 0.0 ? 0 : x > 0.0 ? 1 : -1;
	}

	private double calcError(int index, double a, double k) {
		double s = lnQ[index] + a;
		return W[index] * (1.0 - Math.exp(-s * s / (2.0 * k * k)));
	}

	private double calcErrorDeriv(int index, double a, double k) {
		double s = lnQ[index] + a;
		return W[index] * Math.exp(-s * s / (2.0 * k * k)) * s / (k * k);
	}

	private double calcErrorSum(double a, double k) {
		double sum = 0.0;
		for (int i = 0; i < lnQ.length; i++)
			sum += calcError(i, a, k);

		return sum;
	}

	private double calcErrorDerivSum(double a, double k) {
		double sum = 0.0;

		for (int i = 0; i < lnQ.length; i++)
			sum += calcErrorDeriv(i, a, k);

		return sum;
	}

	private double[][] getMinBounds(double k0, double a0, double d, double s, AceProgress aceProgress, float progressCoeff, float progressCoeffInit) {
		Vector<double[]> v = new Vector<double[]>();

		double[] rl;
		int len = (int) (Math.round(2.0 * s / d));
		int lenP = len / 25;

		double ai;
		int last_sgn, sgn;

		int intersz = 0;

		ai = a0 - s;
		last_sgn = SGN(calcErrorDerivSum(ai, k0));

		for (int i = 1; i <= len; i++) {

			ai = a0 - s + 2.0 * s * i / len;

			sgn = SGN(calcErrorDerivSum(ai, k0));

			if (sgn == 1 && last_sgn == -1) { // zero crossing -+ detected
				rl = new double[2];
				rl[0] = a0 - s + 2.0 * s * (i - intersz - 1) / len;
				rl[1] = ai;

				v.add(rl);

				if (aceProgress != null)
					aceProgress.setStringOv("found minimum " + v.size());
			}

			if (sgn != 0) {
				last_sgn = sgn;
				intersz = 0;
			} else intersz++;

			if (aceProgress != null && lenP != 0 && i % lenP == 0)
				aceProgress.setPartOv(progressCoeffInit + ((double) (i + 1) / (double) len) * ((progressCoeff - progressCoeffInit) / 2.0f));
		}

		// Convert vector into array
		double ret[][] = new double[v.size()][2];
		for (int i = 0; i < v.size(); i++) {
			ret[i][0] = ((v.elementAt(i)))[0];
			ret[i][1] = ((v.elementAt(i)))[1];
		}

		return ret.length == 0 ? null : ret;

	}

	private double minimizeError(double aMin, double aMax, double k0, double convMax) {
		int count = 0;
		int maxIter = 400;
		double k = k0;

		double loA = aMin; // aIn-inter;
		double hiA = aMax; // aIn+inter;

		double loD = calcErrorDerivSum(loA, k);
		double hiD = calcErrorDerivSum(hiA, k);

		double ex, inA, inD;

		if (loD > hiD) {
			ex = loD;
			loD = hiD;
			hiD = ex;
			ex = loA;
			loA = hiA;
			hiA = ex;
		}

		if (loD * hiD > 0.0)
			System.out.println("!!!!!!!!!!!!!!!!!!!");

		do {
			inA = (loA + hiA) * 0.5;
			inD = calcErrorDerivSum(inA, k);

			if (inD > 0.0) {
				hiD = inD;
				hiA = inA;
			} else {
				loD = inD;
				loA = inA;
			}
			count++;
		}
		while (count < maxIter && (inD > convMax || inD < -convMax));
		return inA;
	}

	public void setInput(AceTable x, AceTable y) {
		System.out.println("Normalization of " + name);
		x.sortByProbe();
		y.sortByProbe();

		TableRow xd[] = x.getAsTableRow();
		TableRow yd[] = y.getAsTableRow();

		int sizex = x.getNumRows();
		int sizey = y.getNumRows();

		int curx = 0;
		int cury = 0;

		int match = 0;

		while (curx < sizex && cury < sizey) {
			while (xd[curx].probeID.compareTo(yd[cury].probeID) != 0) {
				while (cury < sizey && xd[curx].probeID.compareTo(yd[cury].probeID) < 0)
					curx++;
				if (curx >= sizex)
					break;
				while (cury < sizey && xd[curx].probeID.compareTo(yd[cury].probeID) > 0)
					cury++;
				if (cury >= sizey)
					break;
			}
			if (curx >= sizex || cury >= sizey)
				break;

			if (xd[curx].signal > 0.0f && xd[curx].variation > 0.0f && yd[cury].signal > 0.0f && yd[cury].variation > 0.0f)
				match++;

			curx++;
			cury++;
		}

		lnQ = new float[match];
		W = new float[match];

		match = 0;
		curx = 0;
		cury = 0;

		while (curx < sizex && cury < sizey) {
			while (xd[curx].probeID.compareTo(yd[cury].probeID) != 0) {
				while (curx < sizex && xd[curx].probeID.compareTo(yd[cury].probeID) < 0)
					curx++;
				if (curx >= sizex)
					break;
				while (cury < sizey && xd[curx].probeID.compareTo(yd[cury].probeID) > 0)
					cury++;
				if (cury >= sizey)
					break;
			}
			if (curx >= sizex || cury >= sizey)
				break;

			if (xd[curx].probeID.compareTo(yd[cury].probeID) != 0)
				System.out.println("not matching ids");

			if (xd[curx].signal > 0.0f && xd[curx].variation > 0.0f && yd[cury].signal > 0.0f && yd[cury].variation > 0.0f) {
				lnQ[match] = (float) (Utils.log2(xd[curx].signal / yd[cury].signal));
				W[match] = (float) (1.0 / Math.sqrt(xd[curx].variation * xd[curx].variation + yd[cury].variation * yd[cury].variation));
				match++;
			}

			curx++;
			cury++;
		}
	}

	public void setInput(AceTable x, AceTable y, String name) {
		setName(name);
		setInput(x, y);
	}

	/**
	 * Calculates the NeONORM bias coefficient for the previously specified
	 * AceTables
	 * 
	 * @param aIn
	 *             center of search interval
	 * @param inter
	 *             half width of search interval
	 * @param k0
	 *             k koefficient of NeONORM function
	 * @param convMax
	 *             error differnece threshold
	 * @return the NeoNORM bias coefficient
	 */
	public double absoluteMinimum(double aIn, double inter, double k0, double convMax) {
		return absoluteMinimum(aIn, inter, k0, convMax, null, 1, 1);
	}

	public double absoluteMinimum(double aIn, double inter, double k0, double convMax, AceProgress aceProgress, float progressCoeff, float progressCoeffInit) {
		if (k0 > 2.0)
			return minimizeError(aIn - inter, aIn + inter, k0, convMax);

		double[][] cand = getMinBounds(k0, aIn, k0 / 5.0, inter, aceProgress, progressCoeff, progressCoeffInit);

		if (cand == null)
			return Double.NaN;


		double[] min_a = new double[cand.length];
		double[] min_e = new double[cand.length];
		for (int i = 0; i < cand.length; i++) {
			min_a[i] = minimizeError(cand[i][0], cand[i][1], k0, convMax);
			min_e[i] = calcErrorSum(min_a[i], k0);
			if (aceProgress != null)
				aceProgress.setPartOv((i + 1) / cand.length * progressCoeff);
		}

		double abs_min_a = min_a[0];
		double abs_min_e = min_e[0];

		for (int i = 1; i < min_a.length; i++) {
			if (abs_min_e > min_e[i]) {
				abs_min_e = min_e[i];
				abs_min_a = min_a[i];
			}
		}

	}

	public double minimizeErrorPlot(double aIn, double inter, double k0, double convMax, File f, AceProgress aceProgress, float progressCoeff) {
		int count = 0;
		int maxIter = 400;

		double minE = -3.0;
		double maxE = 2.0;

		double k = k0;

		double loA = aIn - inter;
		double hiA = aIn + inter;

		double loD = calcErrorDerivSum(loA, k);
		double hiD = calcErrorDerivSum(hiA, k);

		double ex;

		if (loD > hiD) {
			ex = loD;
			loD = hiD;
			hiD = ex;
			ex = loA;
			loA = hiA;
			hiA = ex;
		}

		if (loD * hiD > 0.0)
			System.out.println("!!!!!!!!!!!!!!!!!!!");

		double inA, inD;

		FileWriter fw = null;
		if (f != null) {
			try {
				fw = new FileWriter(f);
			}
			catch (Exception e1) {}
		}

		do {
			inA = (loA + hiA) * 0.5;
			inD = calcErrorDerivSum(inA, k);

			if (inD > 0.0) {
				hiD = inD;
				hiA = inA;
			} else {
				loD = inD;
				loA = inA;
			}
			count++;

			if (fw != null) {
				try {
					fw.write(Utils.getDecs(inA, 6) + "\t" + Utils.getDecs(calcErrorSum(inA, k), 6) + Utils.lineEnd);
				}
				catch (Exception e2) {}
			}
		}
		while (count < maxIter && (inD > convMax || inD < -convMax));

		if (fw != null) {
			try {
				fw.close();
			}
			catch (Exception e3) {}
		}

		double tmp, p;

		double rng_a = 4.0;

		int res_a = 200;
		int res_k = 200;

		double dat[][] = new double[res_a + 1][res_k + 1];

		FileWriter fwr;
		try {
			fwr = new FileWriter(AceUtil.getTemp());

			fwr.write("x");

			for (int ik = 0; ik <= res_k; ik++) {
				k = Utils.exp10(minE + (double) ik / (double) res_k * (maxE - minE));
				fwr.write("\t" + Utils.getDecs(k, 4, false));
			}
			fwr.write(Utils.lineEnd);

			for (int i = 0; i < res_a; i++) {
				p = inA + ((((double) i / (double) res_a) - 0.5) * 2.0 * rng_a);

				fwr.write(Utils.getDecs(p, 4, false));

				for (int ik = 0; ik <= res_k; ik++) {
					k = Utils.exp10(minE + (double) ik / (double) res_k * (maxE - minE));

					tmp = calcErrorDerivSum(p, k);
					dat[i][ik] = tmp;

					fwr.write("\t" + Utils.getDecs(tmp, 4, false));
				}

				fwr.write(Utils.lineEnd);

			}

			fwr.close();
		}
		catch (Exception e0) {}

		String comment = new String("a[k] : " + Utils.getDecs(absoluteMinimum(aIn, 3.0, 0.5, 0.0001, null, 1, 1), 3, false) + "[0.5], " + Utils.getDecs(absoluteMinimum(aIn, 3.0, 1.0, 0.0001, null, 1, 1), 3, false) + "[1.0], " + Utils.getDecs(absoluteMinimum(aIn, 3.0, 2.0, 0.0001, null, 1, 1), 3, false) + "[2.0], " + Utils.getDecs(absoluteMinimum(aIn, 3.0, 100.0, 0.0001, null, 1, 1), 3, false) + "[100]     ");

		BMPFrame bmpf = new BMPFrame("k-scan view", (float) (inA - rng_a), (float) (inA + rng_a), (float) minE, (float) maxE, name, comment);

		ImageIcon windowIcon = AceXapConfig.getImage("window_icon.jpg");
		bmpf.setIconImage(windowIcon.getImage());
		bmpf.setSize(580, 580);

		bmpf.setData(dat);

		Utils.mainFrame.toBack();
		bmpf.toFront();

		return inA;
	}
}
