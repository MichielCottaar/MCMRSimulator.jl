module DefineSequence
import ...Scanners: Scanner
import ...Sequences: Sequence
import ..BuildingBlocks: BuildingBlock
const defining_sequence::Ref{Bool} = Ref(false)

function define_sequence(f::Function, scanner::Scanner, TR=nothing)
    if defining_sequence[]
        # scanner has already been set; do nothing
        return f()
    else
        defining_sequence[] = true
        try
            return Sequence(BuildingBlock(f()); TR=TR, scanner=scanner)
        finally
            defining_sequence[] = false
        end
    end
end
end