// In the code that follows, we use Rohan Sawhney's linear algebra library for matrices and complex numbers.
// https://github.com/rohan-sawhney/linear-algebra-js

class HobbyPoint {
    // A class for associating numerical quantities from Hobby's algorithm with points that appear on a Hobby curve.
    constructor(x, y, tension){
        this.x = x;
        this.y = y;
        this.cmplx = new Complex(x, y);
        // In what follows, we use Knuth's notation in our variable names.
        this.alpha = 1 / tension;
        this.beta = 1 / tension;
        this.d_val = 0;  // Distance between this point and next.
        this.theta = 0;  // Angle of polygonal line from this point to next.
        this.phi = 0;  // Offset angle.
        this.psi = 0;  // Another offset angle.
    }
}

class HobbyCurve{
    // A class for calculating the control points required to draw a Hobby curve.
    constructor(points, tension, cyclic, begin_curl, end_curl){
        this.points = [];
        for (const point of points){
            let x = point[0];
            let y = point[1];
            console.log(x,y);
            this.points.push(new HobbyPoint(x,y, tension));
        }
        this.ctrl_pts = [];
        this.is_cyclic = cyclic;
        this.begin_curl = begin_curl;
        this.end_curl = end_curl;
        this.num_points = points.length;
    }

    get_ctrl_points(){
        // Calculates and returns all of the control points of the Hobby curve.
        this.calculate_d_vals();
        this.calculate_psi_vals();
        this.calculate_theta_vals();
        this.calculate_phi_vals();
        const ctrl_pts = this.calculate_ctrl_pts();
        this.ctrl_pts = ctrl_pts;
        return ctrl_pts;
    }

    calculate_d_vals(){
        // Calculates the pairwise distances between the points.
        // Skip last point if path is non-cyclic
        const point_inds = (this.is_cyclic ? range(0, this.num_points) : range(0, this.num_points - 1) );
        for (const i of point_inds){
            let z_i = this.points[i % this.num_points];
            const z_j = this.points[(i + 1) % this.num_points];
            z_i.d_val = (z_i.cmplx.minus(z_j.cmplx)).norm();
            console.log("d_val", z_i.d_val)
        }
    }
    
    calculate_psi_vals(){
        // Calculates the psi values by subtracting pairwise phases.
        // Skip first and last point if path is non-cyclic
        const point_inds = (this.is_cyclic ? range(0, this.num_points) : range(1, this.num_points - 1) );
        for (const i of point_inds){
            let z_h = (i == 0 ? this.points[this.num_points - 1] : this.points[i - 1])
            let z_i = this.points[i];
            let z_j = this.points[(i + 1) % this.points.length];
            const polygonal_turn = (z_j.cmplx.minus(z_i.cmplx)).overComplex(z_i.cmplx.minus(z_h.cmplx));
            z_i.psi = Math.atan2(polygonal_turn.im, polygonal_turn.re);
            console.log("psi", z_i.psi)
        }
    }

    calculate_theta_vals(){
        // Calculates the theta values by creating a linear system whose solutions are the values.
        let A = DenseMatrix.zeros(this.num_points);  // Inappropriate names, but they mirror Knuth's notation.
        let B = DenseMatrix.zeros(this.num_points);
        let C = DenseMatrix.zeros(this.num_points);
        let D = DenseMatrix.zeros(this.num_points);
        let R = DenseMatrix.zeros(this.num_points);

        // Calculate the entries of the five vectors.
        // Skip first and last point if path is non-cyclic
        const point_inds = (this.is_cyclic ? range(0, this.num_points) : range(1, this.num_points - 1) );
        for (const i of point_inds){
            let z_h = ( i == 0 ? this.points[this.num_points - 1] : this.points[i - 1]);
            let z_i = this.points[i];
            let z_j = this.points[(i + 1) % this.num_points];

            A.set(z_h.alpha / (z_i.beta ** 2 * z_h.d_val)     , i);
            B.set((3 - z_h.alpha) / (z_i.beta ** 2 * z_h.d_val), i);
            C.set((3 - z_j.beta) / (z_i.alpha ** 2 * z_i.d_val), i);
            D.set(z_j.beta / (z_i.alpha ** 2 * z_i.d_val)     , i);
            R.set(-B.get(i) * z_i.psi - D.get(i) * z_j.psi, i);


            console.log("A", z_h.alpha / (z_i.beta ** 2 * z_h.d_val)     , i)
            console.log("B", (3 - z_h.alpha) / (z_i.beta ** 2 * z_h.d_val), i)
            console.log("C", (3 - z_j.beta) / (z_i.alpha ** 2 * z_i.d_val), i)
            console.log("D", z_j.beta / (z_i.alpha ** 2 * z_i.d_val)     , i)
            console.log("R", -B.get(i) * z_i.psi - D.get(i) * z_j.psi, i)
        }
        // Special consideration for non-cyclic paths.
        if (!this.is_cyclic){
            A.set(0, 0);
            B.set(0, 0);
            C.set(this.num_points - 1, 0);
            D.set(this.num_points - 1, 0);
        }
        // Set up matrix M such that the soln. Mx = R are the theta values.
        let T = new Triplet(this.num_points, this.num_points);
        for (const i of range(0, this.num_points)){
            // Fill i-th row ofT 
            if (i == 0){
                T.addEntry(A.get(i), i, this.num_points - 1);
                console.log(A.get(i), i, this.num_points - 1);
            }
            else{
                T.addEntry(A.get(i), i, i - 1);
                console.log(A.get(i), i,  i- 1);
            }
            T.addEntry(B.get(i) + C.get(i), i, i);
            T.addEntry(D.get(i), i, (i + 1) % this.num_points);

            console.log(B.get(i) + C.get(i), i, i)
            console.log(D.get(i), i, (i + 1) % this.num_points)
        }

        // Fix first and last rows of M for non-cyclic paths
        if (!this.is_cyclic){
            // First row
            const alpha_0 = this.points[0].alpha;
            const beta_1 = this.points[1].beta;
            const xi_0 = (alpha_0 ** 2 * this.begin_curl) / beta_1 ** 2;
            // We are forced to do this zeroing-out because Javascript linear algebra libraries are horrible,
            // and this library did not supply a method for appending; just +=-ing. 
            // Zero out current values
            T.addEntry(-(B.get(0) + C.get(0)), 0 ,0)
            T.addEntry(-D.get(0), 0, 1)
            // Put correct value
            T.addEntry(alpha_0 * xi_0 + 3 - beta_1, 0, 0);
            T.addEntry((3 - alpha_0) * xi_0 + beta_1, 0, 1);
            R.set(-((3 - alpha_0) * xi_0 + beta_1) * this.points[1].psi, 0);

            console.log(alpha_0 * xi_0 + 3 - beta_1, 0, 0)            
            console.log((3 - alpha_0) * xi_0 + beta_1, 0, 1)
            // Last row
            const alpha_n_1 = this.points[this.num_points-2].alpha;
            const beta_n = this.points[this.num_points-1].beta;
            const xi_n = (beta_n ** 2 * this.end_curl) / alpha_n_1 ** 2;
            // Zero out current values
            T.addEntry(-A.get(this.num_points-1), this.num_points-1, this.num_points-2)
            T.addEntry(-(B.get(this.num_points-1) + C.get(this.num_points-1)), this.num_points-1, this.num_points-1);
            // Put correct value
            T.addEntry((3 - beta_n) * xi_n + alpha_n_1 , this.num_points-1, this.num_points-2);
            T.addEntry((beta_n * xi_n + 3 - alpha_n_1), this.num_points-1, this.num_points -1); 

            console.log((3 - beta_n) * xi_n + alpha_n_1 , this.num_points-1, this.num_points-2)
            console.log((beta_n * xi_n + 3 - alpha_n_1), this.num_points-1, this.num_points -1)
            R.set(0, this.num_points - 1);
        }
        // Solve for theta values.
        let M = SparseMatrix.fromTriplet(T);
        let lu = M.lu();
        const thetas = lu.solveSquare(R);

        for (let i = 0; i < thetas.nRows(); i++){
            console.log("theta", thetas.get(i))
        }
        for (let i = 0; i < R.nRows(); i++){
            console.log("R", R.get(i))
        }

        console.log(M)
        console.log(R)
        for (let i = 0; i < this.num_points; i++){
            this.points[i].theta = thetas.get(i);
        }
    }

    calculate_phi_vals(){
        // Calculates the phi_k values via the relationship theta_k + phi_k _ psi_k = 0.
        for (let point of this.points){
            point.phi = -(point.psi + point.theta);
        }
    }

    calculate_ctrl_pts(){
        // Calculates the Bezier control points from z_i to z_{i+1}.
        let ctrl_pts = [];
        // Skip last point if path is non-cyclic
        const point_inds = (this.is_cyclic ? range(0, this.num_points) : range(0, this.num_points - 1) );
        for (const i of point_inds){
            console.log("Hi")
            const z_i = this.points[i];
            const z_j = this.points[(i + 1) % this.num_points];
            const rho_coeff = (1/3) * z_i.alpha * velocity(z_i.theta, z_j.phi);
            const sigma_coeff = (1/3) * z_j.beta * velocity(z_j.phi, z_i.theta);
            const exponential_a = new Complex(rho_coeff * Math.cos(z_i.theta), rho_coeff * Math.sin(z_i.theta));
            const exponential_b = new Complex(sigma_coeff * Math.cos(-z_j.phi), sigma_coeff * Math.sin(-z_j.phi));

            const ctrl_pt_a = z_i.cmplx.plus(exponential_a.timesComplex(z_j.cmplx.minus(z_i.cmplx)));
            const ctrl_pt_b = z_j.cmplx.minus(exponential_b.timesComplex(z_j.cmplx.minus(z_i.cmplx)));
            ctrl_pts.push([ctrl_pt_a.re, ctrl_pt_a.im]);
            ctrl_pts.push([ctrl_pt_b.re, ctrl_pt_b.im]);
        }
        return ctrl_pts;
    }
}

function range(start, stop){
    // Python-like range function over nonnegative integers
    return [...Array(stop).keys()].slice(start);
}

function velocity(theta, phi){
    // Metafont's velocity function.
    const numerator = 2 + Math.sqrt(2) * (Math.sin(theta) - (1 / 16) * Math.sin(phi)) * (Math.sin(phi) - (1 / 16) * Math.sin(theta)) * (
            Math.cos(theta) - Math.cos(phi));
    const denominator = (1 + (1 / 2) * (Math.sqrt(5) - 1) * Math.cos(theta) + (1 / 2) * (3 - Math.sqrt(5)) * Math.cos(phi));
    return numerator / denominator;
}