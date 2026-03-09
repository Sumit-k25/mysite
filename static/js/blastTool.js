function downloadFile(url, name) {
	var xhr = new XMLHttpRequest();
	xhr.open("GET", url, true);
	xhr.responseType = "blob";

	xhr.onload = function() {
		var urlCreator = window.URL || window.webkitURL;
		var fileUrl = urlCreator.createObjectURL(this.response);
		var tag = document.createElement('a');
		tag.href = fileUrl;
		tag.download = name+'.txt';
		document.body.appendChild(tag);
		tag.click();
		document.body.removeChild(tag);
	};
	xhr.send();
}




var resetBtn = document.querySelector('.btn-outline-primary[type="reset"]');
resetBtn.addEventListener('click', function() {
	var fileInput = document.getElementById('fileInput');
	var label = document.getElementById('filenameLabel');
	fileInput.value = '';
	label.innerHTML = 'No file selected';
});



let isValid=false;
const fileInput = document.getElementById('fileInput');
fileInput.addEventListener('change', (event) => {
	const file = event.target.files[0];
	var label = document.getElementById("filenameLabel");

	if (file.name.split('.').pop() === 'fasta') {
		var filename = file.name;
		label.innerHTML = filename;
		isValid=true;
		
	} else {
		fileInput.value = '';
		label.innerHTML = 'Invalid file format. Upload a .fasta file';
		alert("Invalid file format. Upload a .fasta file");
		isValid=false;
	}
});




function submitForm() {
	if(isValid===true) {
		document.getElementById("dataForm").submit();
		const loaderContainer = document.querySelector('.loading-screen');
		loaderContainer.style.display = 'flex';
	}
}