import org.json.JSONObject;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.math.BigInteger;
import java.util.*;

public class LagrangeInterpolation {

    public static void main(String[] args) throws Exception {
        // Read JSON file
        String content = new String(Files.readAllBytes(Paths.get("input.json")));
        JSONObject json = new JSONObject(content);

        // Extract keys (n, k) - not really used but kept for completeness
        int n = json.getJSONObject("keys").getInt("n");
        int k = json.getJSONObject("keys").getInt("k");

        // Decode points
        List<BigInteger> xValues = new ArrayList<>();
        List<BigInteger> yValues = new ArrayList<>();

        for (String key : json.keySet()) {
            if (!key.equals("keys")) {
                JSONObject obj = json.getJSONObject(key);
                int base = Integer.parseInt(obj.getString("base"));
                String value = obj.getString("value");

                // Decode large numbers safely
                BigInteger decoded = new BigInteger(value, base);

                // Key itself is x (convert from string)
                BigInteger x = new BigInteger(key);

                xValues.add(x);
                yValues.add(decoded);
            }
        }

        // Build polynomial coefficients
        Fraction[] coeffs = buildPolynomial(xValues, yValues);

        // Print polynomial
        System.out.print("f(x) = ");
        for (int i = coeffs.length - 1; i >= 0; i--) {
            Fraction c = coeffs[i];
            if (c.num.equals(BigInteger.ZERO)) continue; // skip zero terms
            if (i != coeffs.length - 1 && c.num.signum() >= 0) System.out.print(" + ");
            if (i != coeffs.length - 1 && c.num.signum() < 0) System.out.print(" - ");

            Fraction absC = new Fraction(c.num.abs(), c.den);

            if (!(absC.den.equals(BigInteger.ONE) && absC.num.equals(BigInteger.ONE) && i != 0)) {
                System.out.print(absC);
            }
            if (i > 0) {
                System.out.print("x");
                if (i > 1) System.out.print("^" + i);
            }
        }
        System.out.println();

        // Constant term f(0)
        System.out.println("Constant c = " + coeffs[0].simplify());
    }

    // Build polynomial coefficients with Lagrange interpolation
    public static Fraction[] buildPolynomial(List<BigInteger> x, List<BigInteger> y) {
        int n = x.size();
        Fraction[] coeffs = new Fraction[n];
        Arrays.fill(coeffs, new Fraction(BigInteger.ZERO, BigInteger.ONE));

        for (int i = 0; i < n; i++) {
            // L_i(x): start with [1]
            Fraction[] basis = {new Fraction(BigInteger.ONE, BigInteger.ONE)};

            Fraction denom = new Fraction(BigInteger.ONE, BigInteger.ONE);

            for (int j = 0; j < n; j++) {
                if (j == i) continue;

                // Multiply basis by (x - xj)
                Fraction[] newBasis = new Fraction[basis.length + 1];
                Arrays.fill(newBasis, new Fraction(BigInteger.ZERO, BigInteger.ONE));

                for (int k = 0; k < basis.length; k++) {
                    // term * x
                    newBasis[k + 1] = newBasis[k + 1].add(basis[k]);
                    // term * (-xj)
                    newBasis[k] = newBasis[k].add(basis[k].multiply(new Fraction(x.get(j).negate(), BigInteger.ONE)));
                }
                basis = newBasis;

                // Denominator *= (xi - xj)
                denom = denom.multiply(new Fraction(x.get(i).subtract(x.get(j)), BigInteger.ONE));
            }

            // Multiply basis by yi/denom
            Fraction scale = new Fraction(y.get(i), BigInteger.ONE).multiply(new Fraction(denom.den, denom.num)); // divide
            for (int k = 0; k < basis.length; k++) {
                coeffs[k] = coeffs[k].add(basis[k].multiply(scale));
            }
        }
        return coeffs;
    }
}

// Helper Fraction class (exact rational arithmetic)
class Fraction {
    BigInteger num, den;

    Fraction(BigInteger num, BigInteger den) {
        if (den.equals(BigInteger.ZERO)) throw new ArithmeticException("Divide by zero");
        if (den.signum() < 0) { // normalize sign
            num = num.negate();
            den = den.negate();
        }
        this.num = num;
        this.den = den;
    }

    Fraction add(Fraction other) {
        BigInteger newNum = this.num.multiply(other.den).add(other.num.multiply(this.den));
        BigInteger newDen = this.den.multiply(other.den);
        return new Fraction(newNum, newDen).simplify();
    }

    Fraction multiply(Fraction other) {
        return new Fraction(this.num.multiply(other.num), this.den.multiply(other.den)).simplify();
    }

    Fraction simplify() {
        BigInteger g = num.gcd(den);
        return new Fraction(num.divide(g), den.divide(g));
    }

    @Override
    public String toString() {
        return den.equals(BigInteger.ONE) ? num.toString() : (num + "/" + den);
    }
}
