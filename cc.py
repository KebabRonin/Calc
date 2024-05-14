from youtube_transcript_api import YouTubeTranscriptApi

# assigning srt variable with the list
# of dictionaries obtained by the get_transcript() function
try:
    srt = YouTubeTranscriptApi.get_transcript("uDeMunAHgL4")
except:
	pass
# prints the result
with open("captions", "wt",encoding="utf-8") as f:
	for s in srt:
		f.write(f"{s['text']}\n")
# print(srt)