require './rs.rb'

def dump(arr)
	puts("\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|c|}")
	puts("\\hline")
	print("\\textbf{rs} & ")
	15.times do |i|
		print("\\textbf{r%01X} & " % i)
	end

	#arr[0] = 0
	#arr[1] = 0

	puts("\\textbf{0F}\\tabularnewline")
	puts("\\hline")
	print("\\textbf{s0} & ")
	i = 0
	row = 1
	arr.each() do |el|
		i += 1
		print("%02X" % el)
		if i == 16
			puts("\\tabularnewline")
			puts("\\hline")

			if row < 0x10
				print("\\textbf{s%01X} & " % row)
			end

			row += 1
			i = 0
		else
			print(" & ")
		end
	end
	puts("\\end{tabular}")
end

gf = GaloisField.new
gflog = gf.gf_log
gfilog = gf.gf_exp
data = [0xbf, 0x75, 0x65, 0x00]
# mul (X) = gfilog(mod(gflog(a)+gflog(b), 0xff))
# add (+) = a xor b
# RS = Ki X data
# Ki = gfilog(0)
# gfilog[gfilog[0] + 
# RS = gfilog(i) X data(i)
# rs1 = gfilog[(gflog[0] + gflog[0xbf]) % 0xff]
rs1 = 0
rs2 = gfilog[gflog[gfilog[0]] + gflog[0x75]]
rs3 = gfilog[gflog[gfilog[1]] + gflog[0x65]]
rs4 = gfilog[gflog[gfilog[2]] + gflog[0x00]]
puts("%02x %02x %02x %02x: %02x" % [rs1, rs2, rs3, rs4, (rs1 ^ rs2 ^ rs3 ^ rs4)])
exit(0)

msgIn = [ 0x40, 0xd2, 0x75, 0x47, 0x76, 0x17, 0x32, 0x06, 0x27, 0x26, 0x96, 0xc6, 0xc6, 0x96, 0x70, 0xec ]

1.upto(16) do |nsymp|
	0.upto(16) do |errors|
		0.upto((nsymp * 2) + msgIn.size - errors) do |offset|
			next if(errors > nsymp)

			nsym = nsymp * 2
			out = gf.rsEncode(msgIn, nsym)

			errors.times() do |i|
				out[i + offset] = 0xa1
			end
			
			#print("offset %3d, nsym %3d, errors %3d: #{dump(out)} - " % [offset, nsym, errors])

			corrected = gf.correctMessage(out, nsym)
			if(corrected == nil)
				puts("Can't fix the message")
				exit(1)
			end

			msgIn.each_index() do |i|
				if(msgIn[i] != corrected[i])
					puts("Nope!")
					exit(1)
				end
			end

			#puts("ok")
		end
	end
end
