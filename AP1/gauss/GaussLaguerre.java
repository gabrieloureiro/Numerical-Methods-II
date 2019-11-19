//#EQUIPE 10
//#GABRIEL FRANÇA LOUREIRO  //#MATRÍCULA: 388835
//#ABNER DE LIMA ARAÚJO    //#MATRÍCULA: 398067
//#MN2 2019
import java.lang.Math;
public class GaussLaguerre {

	public static void main(String[] args) {
		GaussLaguerre gl = new GaussLaguerre();								//Grau
		System.out.println("\nO valor da integral é: " + gl.integralPartitions(3) + "\n"); 
	}

	private double f(double x) {
		// Usando para os testes a função 30x⁴ - 2x³ + 8x - 12.5.
		// A integral dessa função no intervalo de -1 a 1 é -13.
		return Math.sin(x);
	}

	public double laguerrePonto1() {
		double root1 = 0.5857864376269049511983;
		double root2 = 3.414213562373095048802;

		double value = (0.8535533905932737622 * f(root1) + 0.1464466094067262377996 * f(root2));

		return value;
	}

	public double laguerrePonto2() {
		double root1 = 0.4157745567834790833115;
		double root2 = 2.294280360279041719822;
		double root3 = 6.289945082937479196866;

		double value = (0.71109300992917301545 * f(root1) + 0.278517733569240848801 * f(root2)
				+ 0.010389256501586135749 * f(root3));

		return value;
	}

	public double laguerrePonto3() {
		double root1 = 0.3225476896193923118004;
		double root2 = 1.745761101158346575687;
		double root3 = 4.536620296921127983279;
		double root4 = 9.395070912301133129234;

		double value = (0.60315410434163360164 * f(root1) + 0.357418692437799686641 * f(root2)
				+ 0.038887908515005384272 * f(root3) + 5.392947055613274501E-4 * f(root4));

		return value;
	}

	public double laguerrePonto4() {
		double root1 = 0.2635603197181409102031;
		double root2 = 1.413403059106516792218;
		double root3 = 3.596425771040722081223;
		double root4 = 7.085810005858837556922;
		double root5 = 12.64080084427578265943;

		double value = (0.5217556105828086524759 * f(root1) + 0.39866681108317592745 * f(root2)
				+ 0.07594244968170759539 * f(root3) + 0.0036117586799220484545 * f(root4)
				+ 2.3369972385776227891 * Math.pow(10, -5) * f(root5));

		return value;
	}

	public double laguerrePonto5() {
		double root1 = 0.2228466041792606894644;
		double root2 = 1.188932101672623030743;
		double root3 = 2.992736326059314077691;
		double root4 = 5.77514356910451050184;
		double root5 = 9.837467418382589917716;
		double root6 = 15.98287398060170178255;

		double value = (0.45896467394996359357 * f(root1) + 0.417000830772120994113 * f(root2)
				+ 0.11337338207404497574 * f(root3) + 0.010399197453149074899 * f(root4)
				+ 2.6101720281493205948 * Math.pow(10, -4) * f(root5)
				+ 8.9854790642962123883 * Math.pow(10, -7) * f(root6));

		return value;
	}

	public double integralPartitions(int roots) {
		double sum = 0;
		switch (roots) {
		case 2:

			sum = laguerrePonto1();

			break;

		case 3:

			sum = laguerrePonto2();

			break;

		case 4:

			sum = laguerrePonto3();

			break;

		case 5:

			sum = laguerrePonto4();

			break;

		case 6:

			sum = laguerrePonto5();

			break;
		}

		return sum;
	}

}
