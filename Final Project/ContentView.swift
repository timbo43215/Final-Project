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
    @State var N: String = "10.0"
    @State var J: String = "1.0"
    @State var g: String = "1.0"
    @State var B: String = "0.0"
    @State var potentialArray: [Double] = []
    @StateObject var mySpins = Spins()
    @StateObject var myEnergy = Energy()
    @StateObject var myPotential = Potential()
    
    var body: some View {
        VStack {
            TextField("N:", text: $N)
            Button(action: {
                self.calculateColdSpinConfiguration1D()
            }) {
                Text("Calculate Arbitrary Spin Configuration")
            }
            Button(action: {
                self.calculateTrialConfiguration1D()
            }) {
                Text("Calculate Trial Configuration")
            }
            Button(action: {
                self.calculateEnergyOfTrialConfiguration1D()
            }) {
                Text("Calculate Energy of Trial Configuration")
            }
        }
        
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
    
    func calculateColdSpinConfiguration1D () {
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        
        for x in 1...upperLimitInteger {
            let spinValue = upOrDown[0]
            spinArray.append(spinValue)
        }
        print(spinArray)
    }
    
    func calculateTrialConfiguration1D () {
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        let randParticle1D = Int.random(in: 0...(upperLimitInteger - 1))
        if (spinArray[randParticle1D] == 0.5) {
            spinArray[randParticle1D] = -0.5
        }
        else {
            spinArray[randParticle1D] = 0.5
        }
        print(spinArray)
    }
    
    func calculateArbitrarySpinConfiguration1D () {
        let N = Double(N)!
        let upperLimit = pow(2.0, N)
        let upperLimitInteger = Int(upperLimit)
        
        for x in 1...upperLimitInteger {
            let s = Int.random(in: 0...1)
            let spinValue = upOrDown[s]
            spinArray.append(spinValue)
        }
    }
    
    func calculateEnergyOfTrialConfiguration1D () {
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
        let bohrMagneton = (eValue*hbarc)/(2.0*m)
        let B = Double(B)!
        var energy = 0.0
        
        for x in 1...(upperLimitInteger - 1) {
                let energyValue = -J*(spinArray[x-1]*spinArray[x]) - (B*bohrMagneton*spinArray[x])
                energy = energy + energyValue
            }
        print(energy)
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
        finalEnergy = -J*totalTotalEnergy
        print(finalEnergy)
    }
}

// func randomElement() -> Self.Element?

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
