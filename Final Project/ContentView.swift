//
//  ContentView.swift
//  Final Project
//
//  Created by Tim Stack on 4/7/23.
//

import SwiftUI
import Charts

struct ContentView: View {
    
    @State var upOrDown = [0.5, -0.5]
    @State var spinArray: [Double] = []
    @State var nextSpinArray: [Double] = []
    @State var timeArray: [Double] = []
    @State var N: String = "2.0"
    @State var J: String = "1.0"
    @State var g: String = "1.0"
    @State var B: String = "0.0"
    @State var kT: String = "1.0"
    @State var potentialArray: [Double] = []
    @State var trialEnergy: Double = 0.0
    @State var energy: Double = 0.0
    @State var magnetizationForPlot: [Double] = []
    @State var specificHeat: [Double] = []
    @State var kTForX: [Double] = []
    @State var thermalProperties = [ThermalProperties]()
    @StateObject var mySpins = Spins()
    @StateObject var myEnergy = Energy()
    @StateObject var myPotential = Potential()
    @StateObject private var twoDMagnet = TwoDMagnet()
    let upColor = Color(red: 0.25, green: 0.5, blue: 0.75)
    let downColor = Color(red: 0.75, green: 0.5, blue: 0.25)
    
    var body: some View {
        HStack {
            VStack {
                Text("Magnetization vs. kT")
                HStack {
                    Chart {
                        ForEach(thermalProperties, id: \.kT) { item in
                            PointMark(
                                x: .value("kT", item.kT),
                                y: .value("Magnetization", item.magnetization)
                            )
                        }
                    }
                    .padding()
                }
            }
            VStack {
                Text("Specific Heat vs. kT")
                Chart {
                    ForEach(thermalProperties, id: \.kT) { item in
                        PointMark(
                            x: .value("kT", item.kT),
                            y: .value("Specific Heat", item.specificHeat)
                        )
                    }
                }
                .padding()
            }
        }
            HStack {
                VStack {
                    Text("Energy vs. kT")
                    Chart {
                        ForEach(thermalProperties, id: \.kT) { item in
                            PointMark(
                                x: .value("kT", item.kT),
                                y: .value("Energy", item.energy)
                            )
                        }
                    }
                    .padding()
                }
                
                TimelineView(.animation) { timeline in
                    Canvas { context, size in
                        twoDMagnet.update(to: timeline.date, N: Int(Double(N)!), spinConfiguration: mySpins.spinConfiguration, isThereAnythingInMyVariable: false)
                        
                        for spin in twoDMagnet.spins {
                            let N = Double(N)!
                            let upperLimit = pow(2.0, N)
                            let upperLimitInteger = Int(upperLimit)
                            let rect = CGRect(x: spin.x * (size.width/CGFloat(Float(upperLimitInteger))), y: spin.y * (size.height/CGFloat(Float(upperLimitInteger))), width: (size.height/CGFloat(Float(upperLimitInteger))), height: (size.height/CGFloat(Float(upperLimitInteger))))
                            let shape = Rectangle().path(in: rect)
                            if (spin.spin){
                                context.fill(shape, with: .color(upColor))}
                            else{
                                context.fill(shape, with: .color(downColor))
                            }
                        }
                    }
                }
                .background(.black)
                .ignoresSafeArea()
                .padding()
            }
                //                Button(action: {
                //                    self.calculateColdMetropolisAlgorithm1D()
                //                    self.clearParameters ()})
                //                {Text("Calculate Spin Configuration from Cold Initial")}
                //                Button(action: {
                //                    self.calculateArbitraryMetropolisAlgorithm1D()
                //                    self.clearParameters ()})
                //                {Text("Calculate Spin Configuration from Arbitrary Initial")}
                HStack {
                    Text("N:")
                    TextField("N:", text: $N)
                        .padding()
                    Text("kT:")
                    TextField("kT:", text: $kT)
                        .padding()
                    
                    
                    Button("Start from Cold", action: setupSpinsFromCold)
                        .padding()
                    Button("Start from Arbitrary", action: setupSpinsFromArbitrary)
                        .padding()
                }
            }
    
    func setupSpinsFromCold(){
        let N = Double(N)!
        self.clearParameters ()
        self.calculateColdMetropolisAlgorithm1D()
        twoDMagnet.setup(N: Int(N), spinConfiguration: mySpins.spinConfiguration, isThereAnythingInMyVariable: false)
    }
    
    func setupSpinsFromArbitrary(){
        let N = Double(N)!
        self.clearParameters ()
        self.calculateArbitraryMetropolisAlgorithm1D()
        twoDMagnet.setup(N: Int(N), spinConfiguration: mySpins.spinConfiguration, isThereAnythingInMyVariable: false)
    }
    
    func clearParameters () {
        myEnergy.energy1D = []
        mySpins.spinConfiguration = []
        spinArray = []
        nextSpinArray = []
    }
    
    /// 1. Start with an arbitrary spin configuration α(k) = {s1, s2,...,sN }.
    /// This uses loops to create a 2x2 array of +/- 0.5 spin values randomly.
    func calculateArbitrarySpinConfiguration2D () {
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        for y in 1...upperLimitInteger {
            for x in 1...upperLimitInteger {
                let s = Int.random(in: 0...1)
                let spinValue = upOrDown[s]
                spinArray.append(spinValue)
            }
            mySpins.arbitrarySpinConfiguration.append(spinArray)
            mySpins.trialSpinConfiguration.append(spinArray)
            spinArray.removeAll()
        }
        //print(mySpins.arbitrarySpinConfiguration)
    }
    /// 1. Start with an arbitrary spin configuration α(k) = {s1, s2,...,sN }.
    func calculateColdSpinConfiguration1D () {
        mySpins.spinConfiguration = []
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        var spinValue: [Double] = []
        
        for i in 1...(upperLimitInteger) {
            spinValue.append(0.5)
        }
        mySpins.spinConfiguration.append(spinValue)
        print(mySpins.spinConfiguration[0])
    }
    /// 1. Start with an arbitrary spin configuration α(k) = {s1, s2,...,sN }.
    func calculateArbitrarySpinConfiguration1D () {
        mySpins.spinConfiguration = []
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        var spinValue: [Double] = []
        
        for i in 1...upperLimitInteger {
            let s = Int.random(in: 0...1)
            spinValue.append(upOrDown[s])
        }
        mySpins.spinConfiguration.append(spinValue)
        print(mySpins.spinConfiguration[0])
    }
    /// 2. Generate a trial configuration α(k+1) by
    ///     a. picking a particle i randomly and
    ///     b. flipping its spin.1
    /// Takes the arbitrary spin configuration created in calculateArbitrarySpinConfiguration1D or calculateColdSpinConfiguration1D and copies it.
    /// Then it takes the value of a random particle in the spin configuration and flips its spin.
    func calculateTrialSpinConfiguration1D (x: Int) -> [Double] {
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        let randParticle1D = Int.random(in: 0...(upperLimitInteger - 1))
        var trialSpins = mySpins.spinConfiguration[x-1]
        
        if (x > 0) {
            if (trialSpins[randParticle1D] == 0.5) {
                
                trialSpins[randParticle1D] = -0.5
            }
            else {
                trialSpins[randParticle1D] = 0.5
            }
           // print(mySpins.spinConfiguration[x-1])
        }
        return trialSpins
    }
    func calculateEnergyOfTrialConfiguration1D (x: Int, trialSpins: [Double], J: Int) {
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        trialEnergy = 0.0
        
        let J = Double(J)
        let eValue = 2.7182818284590452353602874713
        // hbarc in eV*Angstroms
        let hbarc = 1973.269804
        // mass of electron in eVc^2
        let m = 510998.95000
        let g = Double(g)!
        let bohrMagneton = (eValue*hbarc)/(2.0*m)
        let B = Double(B)!
        if (x > 0) {
            for i in 1...(upperLimitInteger - 1) {
                let trialEnergyValue = -J*(trialSpins[i-1]*trialSpins[i]) - (B*bohrMagneton*trialSpins[i])
                trialEnergy = trialEnergy + trialEnergyValue
            }
        }
    }
    
    func calculateEnergyOfPreviousConfiguration1D (x:Int, J: Int) {
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        energy = 0.0
        
        let J = Double(J)
        let eValue = 2.7182818284590452353602874713
        // hbarc in eV*Angstroms
        let hbarc = 1973.269804
        // mass of electron in eVc^2
        let m = 510998.95000
        // let g = Double(g)!
        let bohrMagneton = (eValue*hbarc)/(2.0*m)
        let B = Double(B)!
        
        if (x > 0) {
            for i in 1...(upperLimitInteger - 1) {
                let energyValue = -J*(mySpins.spinConfiguration[x-1][i]*mySpins.spinConfiguration[x-1][i]) - (B*bohrMagneton*mySpins.spinConfiguration[x-1][i])
                energy = energy + energyValue
            }
        }
       // print(energy)
    }
    func calculateEnergyCheck (x: Int, trialSpins: [Double]) {
            let uniformRandomNumber = Double.random(in: 0...1)
            if (x > 0) {
                if (trialEnergy <= energy) {
                    mySpins.spinConfiguration.append(trialSpins)
                    myEnergy.energy1D.append(trialEnergy)
                    print("Trial Accepted")
                    print(mySpins.spinConfiguration[x])
                }
                else {
                    if (calculateRelativeProbability() >= uniformRandomNumber){
                        mySpins.spinConfiguration.append(trialSpins)
                        myEnergy.energy1D.append(trialEnergy)
                        print("Trial Accepted")
                        print(mySpins.spinConfiguration[x])
                    }
                    else {
                        //mySpins.spinConfiguration.removeLast()
                        mySpins.spinConfiguration.append(mySpins.spinConfiguration.last!)
                        myEnergy.energy1D.append(energy)
                        
                        print("Trial Rejected")
                        print(mySpins.spinConfiguration[x-1])
                    }
                }
            }
        }
        /// This calculates the relative probability from Equation 15.13 on page 395 in Landau.
        ///  R = exp(-deltaE/kT), where k is the boltzmann constant, T is temperature in Kelvin, and
        ///  deltaE = E(trial) - E(previous)
    func calculateRelativeProbability () -> Double {
            
            let deltaE = trialEnergy - energy
            // units =  m2 kg s-2 K-1
            var R = 0.0
            let kTDouble = Double(kT)!
            R = exp(-deltaE/(kTDouble))
          //print(R)
            
            return R
        }
        /// U
    // From Equation 15.15
    // U(T) = <E>
    //
    
    func calculateInternalEnergy2D () -> Double {
        let energyCount = myEnergy.energy1D.count
        var totalEnergy: Double = 0.0
        var internalEnergy: Double = 0.0
        
        for U in 0...(energyCount - 1) {
            totalEnergy = totalEnergy + myEnergy.energy1D[U]
        }
        let energyCountDouble = Double(energyCount)
        internalEnergy = totalEnergy/energyCountDouble
        print("Internal Energy:")
        print(internalEnergy)
        return internalEnergy
    }
    
    // From Equation 15.14
    //       N
    //      ___
    //      \
    // M  =  >  s
    //  j   /    i
    //      ---
    //      i=1
    
    func calculateMagnetization2D () -> Double {
        let N = Double(N)!
        let upperLimit = sqrt(N)
        let upperLimitInteger = Int(upperLimit)
        var magnetization: Double = 0.0
        
        for j in 0...(upperLimitInteger - 1){
            for i in 0...(upperLimitInteger - 1) {
                magnetization = magnetization + mySpins.spinConfiguration[i][j]
            }
        }
        magnetizationForPlot.append(magnetization)
        
        return magnetization
    }
    
    //Equation 15.17:
    //                             _
    //       1                    |     2
    // U  = --- SUM(from t=1 to M)| (E )
    //  2    M                    |_  t
    //
    func calculateU2 () -> Double {
        var U2: Double = 0.0
        var U2Sum: Double = 0.0
        var energyValueForSum: Double = 0.0
        let N = Double(N)!
        let upperLimit = sqrt(N)
        let upperLimitInteger = Int(upperLimit)
        
        //        for i in 1...(upperLimitInteger - 1) {
        //            energyValueForSum = myEnergy.energy1D[i]
        //            U2Sum = pow(energyValueForSum, 2.0)
        //        }
        U2 = U2Sum/N
        
        return U2
    }
    
    func calculateSpecificHeat () -> Double {
        let N = Double(N)!
        let kT = Double(kT)!
        var specificHeatValue: Double = 0.0
        var U2 = calculateU2()
        var U = calculateInternalEnergy2D()
        
        specificHeatValue = (U2 - pow(U, 2.0))/(pow(N, 2.0)*pow(kT, 2.0))
        specificHeat.append(specificHeatValue)
        
        return specificHeatValue
    }
    
        // 15.4.1 Metropolis Algorithm Implementation
    func calculateColdMetropolisAlgorithm1D () {
            let J: Int = 1
            let N = Double(N)!
            let upperLimit = pow(2.0, N)
            let upperLimitInteger = Int(upperLimit)
            var count: [Double] = []
            var timeValue = 0.0
            calculateColdSpinConfiguration1D ()
            kTForX.append(0.0)
          //  for y in 0...(upperLimitInteger - 1) {
                for x in 1...(upperLimitInteger) {
                    var trialSpins = calculateTrialSpinConfiguration1D(x: x)
                    calculateEnergyOfTrialConfiguration1D(x: x, trialSpins: trialSpins, J: J)
                    calculateEnergyOfPreviousConfiguration1D(x: x, J: J)
                    calculateEnergyCheck(x: x, trialSpins: trialSpins)
                    //trialEnergy = 0.0
                    //energy = 0.0
                    let X = Double(x)
                    kTForX.append(X)
                    calculateU2()
                    calculateSpecificHeat()
                    calculateMagnetization2D()
                    calculateInternalEnergy2D()
                }
        
            for i in 0...(myEnergy.energy1D.count - 1) {
                let kT = Double(kT)!
                
                thermalProperties.append(ThermalProperties(kT: kTForX[i], specificHeat: specificHeat[i], magnetization: magnetizationForPlot[i], energy: myEnergy.energy1D[i]))
            }
        
            print(mySpins.spinConfiguration.count)
            print(mySpins.spinConfiguration)
            print(mySpins.timeComponent)
            print(myEnergy.energy1D)
        
        }
    func calculateArbitraryMetropolisAlgorithm1D () {
            let J: Int = 1
            let N = Double(N)!
            let upperLimit = pow(2.0, N)
            let upperLimitInteger = Int(upperLimit)
            var count: [Double] = []
            var timeValue = 0.0
            calculateArbitrarySpinConfiguration1D ()
                for x in 1...(upperLimitInteger) {
                    var trialSpins = calculateTrialSpinConfiguration1D(x: x)
                    calculateEnergyOfTrialConfiguration1D(x: x, trialSpins: trialSpins, J: J)
                    calculateEnergyOfPreviousConfiguration1D(x: x, J: J)
                    calculateEnergyCheck(x: x, trialSpins: trialSpins)
                    //trialEnergy = 0.0
                    //energy = 0.0
                    let X = Double(x)
                    kTForX.append(X)
                    calculateU2()
                    calculateSpecificHeat()
                    calculateMagnetization2D()
                    calculateInternalEnergy2D()
                }
              //  count.append(timeValue)
              //  mySpins.timeComponent.append(count)
             //   mySpins.spinConfiguration.append(mySpins.spinConfiguration[y])
        for i in 0...(myEnergy.energy1D.count - 1) {
            let kT = Double(kT)!
            
            thermalProperties.append(ThermalProperties(kT: kTForX[i], specificHeat: specificHeat[i], magnetization: magnetizationForPlot[i], energy: myEnergy.energy1D[i]))
        }
            print(mySpins.spinConfiguration.count)
            print(mySpins.spinConfiguration)
            print(mySpins.timeComponent)
            print(myEnergy.energy1D)
        }
        
        /// 2. Generate a trial configuration α(k+1) by
        ///     a. picking a particle i randomly and
        ///     b. flipping its spin.1
        /// Takes the arbitrary spin configuration created in calculateArbitrarySpinConfiguration and copies it.
        /// Then it takes the value of a random particle in the 2x2 matrix and flips its spin.
        func calculateTrialConfiguration2D () {
            let N = Double(N)!
            let upperLimit = pow(2.0, N)
            let upperLimitInteger = Int(upperLimit)
            let randParticleX = Int.random(in: 0...(upperLimitInteger - 1))
            let randParticleY = Int.random(in: 0...(upperLimitInteger - 1))
            if (mySpins.trialSpinConfiguration[randParticleX][randParticleY] == 0.5) {
                mySpins.trialSpinConfiguration[randParticleX][randParticleY] = -0.5
            }
            else {
                mySpins.trialSpinConfiguration[randParticleX][randParticleY] = 0.5
            }
            //print(mySpins.trialSpinConfiguration)
        }
        
        /// 3. Calculate the energy Eαtr of the trial configuration.
        ///
        ///                                  N-1                  N
        /// E(αk) = < a(k)|sum {V(i)}|a(k)> = -J sum    {s(i)*s(i+1)} - B (mu)  sum {s(i)}
        ///                                  i=1                                 b     i=1
        /// This calculates the energy of the Trial Configuration using Equation 15.4 in Landau and creates
        /// a 2x2 matrix of the energy values.
        func calculateEnergyOfTrialConfiguration2D () {
            let N = Double(N)!
            let upperLimit = pow(2.0, N)
            let upperLimitInteger = Int(upperLimit)
            let J = Double(J)!
            let eValue = 2.7182818284590452353602874713
            // hbarc in eV*Angstroms
            let hbarc = 1973.269804
            // mass of electron in eVc^2
            let m = 510998.95000
            let g = Double(g)!
            let gbohrMagneton = g*((eValue*hbarc)/(2.0*m))
            var energy = 0.0
            var totalTotalEnergy = 0.0
            var finalEnergy = 0.0
            let B = Double(B)!
            
            for y in 0...(upperLimitInteger - 1) {
                for x in 0...(upperLimitInteger - 1) {
                    let potentialValue = (spinArray[x]*spinArray[x+1]) - (gbohrMagneton*(spinArray[x]*B))
                    energy = energy + potentialValue
                }
                totalTotalEnergy = totalTotalEnergy + energy
            }
            finalEnergy = Double(-J)*totalTotalEnergy
            print(finalEnergy)
        }
}

    
struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
