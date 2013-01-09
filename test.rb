require './rs.rb'

def dump(arr)
	arr.dup.map() { |i| "%02x" % i }.join(" ")
end

gf = GaloisField.new
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
			
			print("offset %3d, nsym %3d, errors %3d: #{dump(out)} - " % [offset, nsym, errors])

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

			puts("ok")
		end
	end
end
