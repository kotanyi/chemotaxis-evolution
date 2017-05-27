public class Gaussian {
    private double width, height, centre;

    public Gaussian(double width, double height, double centre) {
        this.width = width;
        this.height = height;
        this.centre = centre;
    }

    public double getFoodAt(double position) {
        double food = height * Math.exp(-(position - centre) * (position - centre) / (2 * width * width));

        if(food < Parameters.gaussianDetectionThreshold)
            food = 0;

        return food;
    }

    public double getHeight() {
        return height;
    }
}
