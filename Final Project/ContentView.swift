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
    @State var T: String = "1.0"
    @State var potentialArray: [Double] = []
    @State var trialEnergy: Double = 0.0
    @State var energy: Double = 0.0
    @StateObject var mySpins = Spins()
    @StateObject var myEnergy = Energy()
    @StateObject var myPotential = Potential()
    
    var body: some View {
        VStack {
            TextField("N:", text: $N)
            //            Button(action: {
            //                self.calculateColdSpinConfiguration1D(x: x)
            //            }) {
            //                Text("Calculate Cold Spin Configuration")
            //            }
            //            Button(action: {
            //                self.calculateTrialSpinConfiguration1D()
            //            }) {
            //                Text("Calculate Trial Configuration")
            //            }
            //            Button(action: {
            //                self.calculateEnergyOfTrialConfiguration1D()
            //            }) {
            //                Text("Calculate Energy of Trial Configuration")
            //            }
            Button(action: {
                self.calculateMetropolisAlgorithm1D()
                self.clearParameters ()
            }) {
                Text("Calculate Spin Configuration")
            }
        }
        
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
        spinArray = []
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        
        for x in 1...upperLimitInteger {
            let s = Int.random(in: 0...1)
            let spinValue = upOrDown[s]
            spinArray.append(spinValue)
        }
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
    func calculateEnergyOfTrialConfiguration1D (x: Int, trialSpins: [Double]) {
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        trialEnergy = 0.0
        
        let J = Double(J)!
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
    
    func calculateEnergyOfPreviousConfiguration1D (x:Int) {
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        energy = 0.0
        
        let J = Double(J)!
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
                let energyValue = -J*(mySpins.spinConfiguration[x-1][i]*mySpins.spinConfiguration[x-1][i]) - (B*bohrMagneton*mySpins.spinConfiguration[x-1][i])
                energy = energy + energyValue
            }
        }
        print(energy)
    }
        func calculateEnergyCheck (x: Int, trialSpins: [Double]) {
            let uniformRandomNumber = Double.random(in: 0...1)
            if (x > 0) {
                if (trialEnergy <= energy) {
                    mySpins.spinConfiguration.append(trialSpins)
                    myEnergy.energy1D.append(trialEnergy)
                    print("Trial Accepted")
                    print(mySpins.spinConfiguration[x-1])
                }
                else {
                    if (calculateRelativeProbability() >= uniformRandomNumber){
                        mySpins.spinConfiguration.append(trialSpins)
                        myEnergy.energy1D.append(trialEnergy)
                        print("Trial Accepted")
                        print(mySpins.spinConfiguration[x-1])
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
            else {}
        }
        /// This calculates the relative probability from Equation 15.13 on page 395 in Landau.
        ///  R = exp(-deltaE/kT), where k is the boltzmann constant, T is temperature in Kelvin, and
        ///  deltaE = E(trial) - E(previous)
        func calculateRelativeProbability () -> Double {
            
            let deltaE = trialEnergy - energy
            // units =  m2 kg s-2 K-1
            let k = 1.380649*pow(10.0, -23.0)
            let T = Double(T)!
            let kT = 1.0
            var R = 0.0
            R = exp(-deltaE/(kT))
          //print(R)
            
            return R
        }
        /// U
        func calculateInternalEnergy1D () {
            
        }
        
        func calculateMetropolisAlgorithm1D () {
            let N = Double(N)!
            let upperLimit = pow(2.0, N)
            let upperLimitInteger = Int(upperLimit)
            var count: [Double] = []
            var timeValue = 0.0
            calculateColdSpinConfiguration1D ()
            for y in 0...(upperLimitInteger) {
                for x in 1...(upperLimitInteger) {
                    var trialSpins = calculateTrialSpinConfiguration1D(x: x)
                    calculateEnergyOfTrialConfiguration1D(x: x, trialSpins: trialSpins)
                    calculateEnergyOfPreviousConfiguration1D(x: x)
                    calculateEnergyCheck(x: x, trialSpins: trialSpins)
                    //trialEnergy = 0.0
                    //energy = 0.0
                }
                timeValue = Double(y-1) + 1.0
              //  count.append(timeValue)
              //  mySpins.timeComponent.append(count)
             //   mySpins.spinConfiguration.append(mySpins.spinConfiguration[y])
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
    
//    func setObjectWillChange(theObject:PlotClass){
//
//        theObject.objectWillChange.send()
//
//    }
    
    //    func plot(selectedGraphIndex: Int) {
    //
    //        setObjectWillChange(theObject: self.plotData)
    //        myCalculatePlotData.plotDataModel = self.plotData.plotArray[0]
    //
    //        switch  selectedGraphIndex {
    //
    //        case 0:
    //
    //            myCalculatePlotData.theText = "Potential\n"
    //
    //            myCalculatePlotData.setThePlotParameters(color: "Blue", xLabel: "x", yLabel: "Potential", title: "Potential")
    //
    //            //          myCalculatePlotData.resetCalculatedTextOnMainThread()
    //
    //
    //            var thePlotData :[(x: Double, y: Double)] =  []
    //
    //            for i in 0..<myPotentials.x.count {
    //                let X = myPotentials.x[i]
    //                var Y = myPotentials.Potential[i]
    //                if Y > 100 {
    //                    Y = 100
    //                }
    //                thePlotData.append((x: X, y: Y))
    //            }
    //            myCalculatePlotData.appendDataToPlot(plotData: thePlotData)
    //
    //            myCalculatePlotData.plotDataModel!.changingPlotParameters.xMax = myPotentials.x.max()! + 2.0
    //
    //            myCalculatePlotData.plotDataModel!.changingPlotParameters.xMin = myPotentials.x.min()! - 2.0
    //
    //            myCalculatePlotData.plotDataModel!.changingPlotParameters.yMax = 55.0
    //
    //            myCalculatePlotData.plotDataModel!.changingPlotParameters.yMin = myPotentials.Potential.min()! - 2.0
    //
    //            setObjectWillChange(theObject: self.plotData)
    //
    //        case 1:
    //
    //            myCalculatePlotData.theText = "Functional\n"
    //
    //            myCalculatePlotData.setThePlotParameters(color: "Blue", xLabel: "x", yLabel: "Functional", title: "Fuctional")
    //
    //            //         myCalculatePlotData.resetCalculatedTextOnMainThread()
    //
    //
    //            var thePlotData :[(x: Double, y: Double)] =  []
    //
    //            for i in 0..<myFunctionals.Energy.count {
    //                let X = myFunctionals.Energy[i]
    //                var Y = myFunctionals.Functional[i]
    //                if Y > 100 {
    //                    Y = 100
    //                }
    //                else if Y < -100 {
    //
    //                    Y = -100
    //
    //                }
    //                thePlotData.append((x: X, y: Y))
    //            }
    //
    //
    //        }
    //    }
    // func randomElement() -> Self.Element?
struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
