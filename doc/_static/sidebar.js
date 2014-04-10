$(window).scroll(function() {
    var st = $(".content-wrapper").position().top;
    var wst = $(window).scrollTop();
    $(".sidebar").css({position: 'relative', top: Math.max(0, wst - st)+'px' });
});
