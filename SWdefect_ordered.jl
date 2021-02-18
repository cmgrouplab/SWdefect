import Statistics
using Random, Distributions
using Plots
using DelimitedFiles 

const NUMCELL = 6
LatticeVector = 5.8277997971
pValue = 0.05  #floor(pValue*2*NUMCELL*NUMCELL)
const PI = 3.1415926
T = 1e-6
uniform = Uniform(0, 1)

mutable struct Cell
    atomsPt::Array
    atomsP::Array
end

mutable struct AtomPt
    position :: Array{Float64,1}
    connectivity :: Array
    index :: Int64
end

mutable struct AtomP
    position :: Array{Float64,1}
    connectivity :: Array
    ifRotate :: Bool
    index :: Int64
end 

function newCell(i::Int, j::Int)
    atomPt1 = AtomPt([(i-1)*LatticeVector,(j-1)*LatticeVector],[],2*(j+(i-1)*NUMCELL-1)+1)
    atomPt2 = AtomPt([(0.5+i-1)*LatticeVector,(0.5+j-1)*LatticeVector],[],2*(j+(i-1)*NUMCELL-1)+2)
    atomP1 = AtomP([(0.374085009+i-1)*LatticeVector,(0.125915006+j-1)*LatticeVector],[],false,4*(j+(i-1)*NUMCELL-1)+1)
    atomP2 = AtomP([(0.125915006+i-1)*LatticeVector,(0.625914991+j-1)*LatticeVector],[],false,4*(j+(i-1)*NUMCELL-1)+2)
    atomP3 = AtomP([(0.8740850096+i-1)*LatticeVector,(0.374085009+j-1)*LatticeVector],[],false,4*(j+(i-1)*NUMCELL-1)+3)
    atomP4 = AtomP([(0.625914991+i-1)*LatticeVector,(0.874085009+j-1)*LatticeVector],[],false,4*(j+(i-1)*NUMCELL-1)+4)

    return Cell([atomPt1, atomPt2], [atomP1, atomP2, atomP3, atomP4])
end

function boundaryDisplacement(a, b)
    delta = a - b
    if abs(delta) > LatticeVector
        if delta > 0
            delta -= NUMCELL*LatticeVector
        else
            delta += NUMCELL*LatticeVector
        end
    end
    return delta
end

function boundaryDistance(a, b)
    deltax = boundaryDisplacement(a.position[1], b.position[1])
    deltay = boundaryDisplacement(a.position[2], b.position[2])

    return sqrt(deltax * deltax + deltay * deltay)
end

function initializationCells()
    cells = Array{Cell}(undef, NUMCELL,NUMCELL)
    for i = 1:NUMCELL
        for j = 1:NUMCELL
            cells[i,j] = newCell(i, j)
        end
    end

    for i = 1:NUMCELL
        for j= 1:NUMCELL
            iminus1 = i-1
            jminus1 = j-1
            iplus1 = i+1
            jplus1 = j+1
            if iminus1 == 0
                iminus1 = NUMCELL
            end

            if jminus1 == 0
                jminus1 = NUMCELL
            end

            if iplus1 == NUMCELL + 1
                iplus1 = 1
            end

            if jplus1 == NUMCELL + 1
                jplus1 = 1
            end

            #cells[i,j].atomsPt[1].connectivity = [cells[i,j].atomsP[1],cells[iminus1,j].atomsP[3],cells[i,jminus1].atomsP[2],cells[iminus1,jminus1].atomsP[4]]
            #cells[i,j].atomsPt[2].connectivity = [cells[i,j].atomsP[1],cells[i,j].atomsP[2],cells[i,j].atomsP[3],cells[i,j].atomsP[4]]
            cells[i,j].atomsP[1].connectivity = [cells[i,jminus1].atomsP[4], cells[i,j].atomsPt[1], cells[i,j].atomsPt[2]]
            cells[i,j].atomsP[2].connectivity = [cells[iminus1,j].atomsP[3],cells[i,j].atomsPt[2],cells[i,jplus1].atomsPt[1]]
            cells[i,j].atomsP[3].connectivity = [cells[iplus1,j].atomsP[2],cells[i,j].atomsPt[2],cells[iplus1,j].atomsPt[1]]
            cells[i,j].atomsP[4].connectivity = [cells[i,jplus1].atomsP[1],cells[i,j].atomsPt[2],cells[iplus1,jplus1].atomsPt[1]]
        end
    end
    return cells
end



function rotateAtoB(vectorA::Array, rotateAngle::Number)
    absVectorA = sqrt(vectorA[1] * vectorA[1] + vectorA[2] * vectorA[2])
    vectorAangle = atan(vectorA[2], vectorA[1])
    vectorAangle += rotateAngle 
    vectorA = [absVectorA * cos(vectorAangle),absVectorA * sin(vectorAangle)]
    return vectorA
end

function rotate90(a,b)
    sentry = [0,0,0,0]
    if abs(a[1]-b[1])>LatticeVector
        if a[1] - b[1] < 0
            b[1] = b[1] - NUMCELL * LatticeVector
            sentry[1] = 1
        else
            a[1] = a[1] - NUMCELL * LatticeVector
            sentry[3] = 1
        end
    end

    if abs(a[2]-b[2])>LatticeVector
        if a[2] - b[2] < 0
            b[2] = b[2] - NUMCELL * LatticeVector
            sentry[2] = 1
        else
            a[2] = a[2] - NUMCELL * LatticeVector
            sentry[4] = 1
        end
    end
         

    centerVector = (a+b)/2
    rotateVectorA = a - centerVector
    rotateVectorA = rotateAtoB(rotateVectorA,-PI/2)
    a = centerVector + rotateVectorA

    rotateVectorB = b - centerVector
    rotateVectorB = rotateAtoB(rotateVectorB,-PI/2)
    b = centerVector + rotateVectorB

    if sentry[1] == 1
        b[1] = b[1] +  NUMCELL * LatticeVector
    end

    if sentry[2] == 1
       b[2] = b[2] +  NUMCELL * LatticeVector
    end
    if sentry[3] == 1
        a[1] = a[1] +  NUMCELL * LatticeVector
    end
    if sentry[4] == 1
        a[2] = a[2] +  NUMCELL * LatticeVector
    end
    return a,b
end

function changeConnectivityK1(cells,i,j)
    iminus1 = i-1
    jminus1 = j-1
    iplus1 = i+1
    jplus1 = j+1
    if iminus1 == 0
        iminus1 = NUMCELL
    end

    if jminus1 == 0
        jminus1 = NUMCELL
    end

    if iplus1 == NUMCELL + 1
        iplus1 = 1
    end

    if jplus1 == NUMCELL + 1
        jplus1 = 1
    end
    for index in 2:3
        if cells[i,j].atomsP[1].connectivity[index] == cells[i,j].atomsPt[1]
            cells[i,j].atomsP[1].connectivity[index] = cells[iplus1,j].atomsPt[1]
        end
        
        if  cells[i,j].atomsP[1].connectivity[1].connectivity[index] == cells[iplus1,j].atomsPt[1]
            cells[i,j].atomsP[1].connectivity[1].connectivity[index] = cells[i,j].atomsPt[1]

        end
    end
end

function changeConnectivityK2(cells,i,j)    
    iminus1 = i-1
    jminus1 = j-1
    iplus1 = i+1
    jplus1 = j+1
    if iminus1 == 0
        iminus1 = NUMCELL
    end

    if jminus1 == 0
        jminus1 = NUMCELL
    end

    if iplus1 == NUMCELL + 1
        iplus1 = 1
    end

    if jplus1 == NUMCELL + 1
        jplus1 = 1
    end
    for index in 2:3
        if cells[i,j].atomsP[2].connectivity[index] == cells[i,jplus1].atomsPt[1]
            cells[i,j].atomsP[2].connectivity[index] = cells[i,j].atomsPt[1]
        end
        
        if  cells[i,j].atomsP[2].connectivity[1].connectivity[index] == cells[i,j].atomsPt[1]
            cells[i,j].atomsP[2].connectivity[1].connectivity[index] = cells[i,jplus1].atomsPt[1]

        end
    end
end

function SWdefect(cells)
    for n =1:round(pValue*NUMCELL*NUMCELL)
        println("n: ",n)
        i = rand(1:NUMCELL)
        j = rand(1:NUMCELL)
        k = rand(2:3)

        while cells[i,j].atomsP[k].ifRotate 
            i = rand(1:NUMCELL)
            j = rand(1:NUMCELL)
            k = rand(2:3)
        end

        iminus1 = i-1
        jminus1 = j-1
        iplus1 = i+1
        jplus1 = j+1
        if iminus1 == 0
            iminus1 = NUMCELL
        end
    
        if jminus1 == 0
            jminus1 = NUMCELL
        end
    
        if iplus1 == NUMCELL + 1
            iplus1 = 1
        end
    
        if jplus1 == NUMCELL + 1
            jplus1 = 1
        end

        cells[i,j].atomsP[k].connectivity[1].position,cells[i,j].atomsP[k].position = rotate90(cells[i,j].atomsP[k].connectivity[1].position,cells[i,j].atomsP[k].position)
        cells[i,j].atomsP[k].connectivity[1].ifRotate = true
        cells[i,j].atomsP[k].ifRotate = true
        # if k == 1 
        #     changeConnectivityK1(cells,i,j)
        # end
        
        # if k ==4
        #     changeConnectivityK1(cells,i,jplus1)
        # end

        if k==2
            changeConnectivityK2(cells,i,j)
        end

        if k==3
            changeConnectivityK2(cells,iplus1,j)
        end         

    end
end

function getEnergy(cells)
    PPdistance = []
    PtPdistance = []
    energy1 = 0
    energy2 = 0
    for i = 1:NUMCELL
        for j = 1:NUMCELL
            append!(PPdistance,boundaryDistance(cells[i,j].atomsP[1],cells[i,j].atomsP[1].connectivity[1]))
            append!(PPdistance,boundaryDistance(cells[i,j].atomsP[2],cells[i,j].atomsP[2].connectivity[1]))
            append!(PtPdistance,boundaryDistance(cells[i,j].atomsP[1],cells[i,j].atomsP[1].connectivity[2]))
            append!(PtPdistance,boundaryDistance(cells[i,j].atomsP[1],cells[i,j].atomsP[1].connectivity[3]))
            append!(PtPdistance,boundaryDistance(cells[i,j].atomsP[2],cells[i,j].atomsP[2].connectivity[2]))
            append!(PtPdistance,boundaryDistance(cells[i,j].atomsP[2],cells[i,j].atomsP[2].connectivity[3]))
            append!(PtPdistance,boundaryDistance(cells[i,j].atomsP[3],cells[i,j].atomsP[3].connectivity[2]))
            append!(PtPdistance,boundaryDistance(cells[i,j].atomsP[3],cells[i,j].atomsP[3].connectivity[3]))
            append!(PtPdistance,boundaryDistance(cells[i,j].atomsP[4],cells[i,j].atomsP[4].connectivity[2]))
            append!(PtPdistance,boundaryDistance(cells[i,j].atomsP[4],cells[i,j].atomsP[4].connectivity[3]))
        end
    end
    energy1 = Statistics.std(PPdistance)
    energy2 = Statistics.std(PtPdistance)
    println(PtPdistance)
    return energy1+energy2
end

function evolvePt(cells,i,j,k,x,y)
    cells[i,j].atomsPt[k].position += [x,y]
end

function evolveP(cells,i,j,k,x,y)
    cells[i,j].atomsP[k].position += [x,y]
end

function isAccepted(energy,newEnergy)
    if newEnergy > energy
        deltaEnergy = newEnergy - energy 
        #println("newEnergy: ", newEnergy)
        #println("energy: ",energy) 
        #println("deltaEnergy: ", deltaEnergy) 
        relativePossibility = exp(-deltaEnergy/T)  
        #println("relativePossibility: ", relativePossibility)
        test = rand(uniform)
        #println("test: ",test)
	    if test > relativePossibility
            return false
        end
    end
    return true
end

function relax(cells,numOfRun)
    for run = 1:numOfRun
        energy = getEnergy(cells)
        if run ==1
            println("initial energy:", energy*1e8)
        end

        i = rand(1:NUMCELL)
        j = rand(1:NUMCELL)
        k = rand(1:6)
        x = 0.0001*LatticeVector*(rand(uniform)-0.5)
        y = 0.0001*LatticeVector*(rand(uniform)-0.5)
        if k <3
            evolvePt(cells,i,j,k,x,y)
        else
            evolveP(cells,i,j,k-2,x,y)
        end

        newEnergy = getEnergy(cells)
        if !(isAccepted(energy,newEnergy))
            if k<3
                evolvePt(cells,i,j,k,-x,-y)
            else
                evolveP(cells,i,j,k-2,-x,-y)
            end
        end
        if run == numOfRun
            println("final energy: ", getEnergy(cells)*1e8 )
        end   
    end
end     

function drawPng(cells,pValue,simID)
    xPt = []
    yPt = []
    xP =[]
    yP = []
    for i = 1:NUMCELL
        for j=1:NUMCELL
            append!(xPt,[cells[i,j].atomsPt[1].position[1],cells[i,j].atomsPt[2].position[1]])
            append!(yPt,[cells[i,j].atomsPt[1].position[2],cells[i,j].atomsPt[2].position[2]])
        end
    end

    for i = 1:NUMCELL
        for j=1:NUMCELL
            append!(xP,[cells[i,j].atomsP[1].position[1],cells[i,j].atomsP[2].position[1],cells[i,j].atomsP[3].position[1],cells[i,j].atomsP[4].position[1]])
            append!(yP,[cells[i,j].atomsP[1].position[2],cells[i,j].atomsP[2].position[2],cells[i,j].atomsP[3].position[2],cells[i,j].atomsP[4].position[2]])
        end
    end

    plot(legend=false, aspect_ratio=:equal,grid = false,ticks = false)
    scatter!(xPt, yPt, markersize=8,color=:black)
    scatter!(xP, yP, markersize=6,color=:red)

    for i =1:NUMCELL
        for j=1:NUMCELL
            for aPt in cells[i, j].atomsPt
                for a in aPt.connectivity
                    if sqrt((aPt.position[1]-a.position[1])^2 + (aPt.position[2]-a.position[2])^2) <5
                        plot!([aPt.position[1], a.position[1]], [aPt.position[2], a.position[2]], color=:black)
                    end
                end
            end

            for aP in cells[i, j].atomsP
                for a in aP.connectivity
                    if sqrt((aP.position[1]-a.position[1])^2 + (aP.position[2]-a.position[2])^2) <5
                        if typeof(a) == AtomP
                            if mod(a.index,4) == 1 || mod(a.index,4) == 0
                                plot!([aP.position[1], a.position[1]], [aP.position[2], a.position[2]], color=:red,linewidth= 8)
                            else
                                plot!([aP.position[1], a.position[1]], [aP.position[2], a.position[2]], color=:blue,linewidth= 8)
                            end

                        else
                            plot!([aP.position[1], a.position[1]], [aP.position[2], a.position[2]], color=:black)
                        end
                    end
                end
            end
        end
    end
    savefig("p_$(pValue)/$(simID)/ordered_p_$(pValue)_$(simID).png")
end

function simulation(pValue, simID)
    mkdir("p_$(pValue)/$(simID)")
    pwd()
    cells = initializationCells()
    SWdefect(cells)
    #relax(cells,10)
    positiondata = Array{Union{Int, Float64}}(undef, cells[end,end].atomsPt[end].index, 3)
    positiondata2 = Array{Union{Int, Float64}}(undef, cells[end,end].atomsP[end].index, 6)
    for i = 1:NUMCELL, j = 1:NUMCELL
        for aPt in cells[i, j].atomsPt
            positiondata[aPt.index, 1] = aPt.index
            positiondata[aPt.index, 2:end] = aPt.position
        end
    end
    writedlm("p_$(pValue)/$(simID)/positionPt.csv", positiondata, ',')

    for i = 1:NUMCELL, j=1:NUMCELL
        for aP in cells[i, j].atomsP
            positiondata2[aP.index, 1] = aP.index
            positiondata2[aP.index, 2:3] = aP.position
            positiondata2[aP.index, 4:end] = [con.index for con in aP.connectivity]
        end
    end
    writedlm("p_$(pValue)/$(simID)/positionP.csv", positiondata2, ',')
    drawPng(cells,pValue,simID)   
end

function runSimulation()
    while pValue < 1.01
        mkdir("p_$(pValue)")
        for simID = 6:6
            simulation(pValue, simID)
        end
        global pValue += 0.05
        global pValue = round(pValue, digits=2)
    end
end

@time runSimulation()






