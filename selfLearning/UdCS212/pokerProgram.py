##------------------------------write a poker program

def poker(hands):
	"Return the best hand: poker([hand,...]) => hand"
	return max(hands,key=hand_rank)
print max([3,4,5,0], max([3,4,-5,0], key=abs))

def hand_rank(hand):
	"Return a value indicating the ranking of a hand."
	ranks = card_ranks(hand)
	if straight(ranks) and flush(hand):
		return (8, max(ranks))
	elif kind(4, ranks):
		return (7, kind(4,ranks),kind(1,ranks)) 
	elif ...	

def test():
	"Test cases for the funciton in poker program."
	sf = "6C 7C 8C 9C TC".split()
	fk = "6C 7C 8C 9C TC".split()
	fh = "6C 7C 8C 9C TC".split()
